"""Run-owned, immutable JSON checkpoints with atomic publication and a process lock."""
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import tempfile


def _fsync_directory(path):
    """Flush a directory entry so a linked file survives a crash."""
    directory = os.open(path, os.O_RDONLY)
    try:
        os.fsync(directory)
    finally:
        os.close(directory)


class Checkpoint:
    """Retain completed stages; an explicit resume must match execution identity.

    Incomplete temporary files are preserved and ignored. The OS releases flock
    after interruption, so a crash never requires deleting a stale lock file.
    """

    def __init__(self, path, identity, *, resume=False):
        """Bind a checkpoint directory to an execution identity; resume requires an identical manifest."""
        # Staged files come from tempfile.mkdtemp, which returns absolute paths;
        # a cwd-relative prefix must still map onto the same directory.
        self.path = Path(path).absolute()
        self.identity = identity
        self.resume = resume
        self.lock = None

    def __enter__(self):
        if self.path.is_symlink():
            raise ValueError('Checkpoint directory must not be a symlink')
        if not self.resume:
            self.path.mkdir(parents=True, exist_ok=False)
        elif not self.path.is_dir():
            raise ValueError(f'No checkpoint to resume: {self.path}')
        try:
            lock_path = self.path / '.lock'
            fd = os.open(lock_path, os.O_CREAT | os.O_RDWR | os.O_NOFOLLOW, 0o600)
            self.lock = os.fdopen(fd, 'a')
            fcntl.flock(self.lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
            expected = {'schema_version': 1, 'identity': self.identity}
            if self.resume:
                if self.read('manifest') != expected:
                    raise ValueError('Checkpoint identity mismatch: config, inputs, code or environment changed')
            else:
                self.save('manifest', expected)
            return self
        except BaseException:
            self.__exit__(None, None, None)
            raise

    def __exit__(self, *_):
        if self.lock is not None:
            self.lock.close()
            self.lock = None

    def _file(self, key):
        if not re.fullmatch(r'[A-Za-z0-9_.-]+', key) or key in ('.', '..'):
            raise ValueError('Invalid checkpoint record name')
        path = self.path / (key + '.json')
        if path.is_symlink():
            raise ValueError('Checkpoint records must not be symlinks')
        return path

    def read(self, key):
        """Return a completed record (checksum verified) or None when it does not exist."""
        path = self._file(key)
        if not path.exists():
            return None
        try:
            envelope = json.loads(path.read_text(), parse_constant=self._invalid)
            if not isinstance(envelope, dict) or set(envelope) != {'sha256', 'value'}:
                raise ValueError(f'Malformed checkpoint envelope: {path}')
            encoded = json.dumps(envelope['value'], sort_keys=True, allow_nan=False).encode()
            if hashlib.sha256(encoded).hexdigest() != envelope['sha256']:
                raise ValueError(f'Checkpoint record checksum mismatch: {path}')
            return envelope['value']
        except (json.JSONDecodeError, UnicodeError) as exc:
            raise ValueError(f'Malformed checkpoint record: {path}') from exc

    @staticmethod
    def _invalid(value):
        raise ValueError(f'Nonfinite checkpoint value: {value}')

    def save(self, key, value):
        """Persist a completed record exactly once with fsync and an exclusive hard link."""
        if self.lock is None:
            raise RuntimeError('Checkpoint writes require an active context')
        path = self._file(key)
        canonical = json.dumps(value, sort_keys=True, allow_nan=False).encode()
        envelope = {'sha256': hashlib.sha256(canonical).hexdigest(), 'value': value}
        encoded = json.dumps(envelope, indent=2, sort_keys=True, allow_nan=False) + '\n'
        if path.exists():
            if self.read(key) == value:
                return
            raise FileExistsError(f'Completed checkpoint record already exists: {path}')
        fd, temporary = tempfile.mkstemp(prefix=f'.{key}.', suffix='.pending', dir=self.path)
        with os.fdopen(fd, 'w', encoding='utf-8') as stream:
            stream.write(encoded)
            stream.flush()
            os.fsync(stream.fileno())
        # Exclusive link prevents replacement, including a file appearing mid-write.
        os.link(temporary, path)
        # Keep the completed temporary inode as an artifact; no destructive cleanup.
        _fsync_directory(self.path)


def execution_identity(metadata, *, no_pdf):
    """Identity fields of a run: config, inputs, waveforms, source hashes, packages, platform, threads, PDF mode."""
    keys = ('config_sha256', 'datasets', 'waveforms', 'source_sha256',
            'packages', 'python', 'platform', 'threads')
    return {**{key: metadata[key] for key in keys}, 'no_pdf': bool(no_pdf)}


def publish_outputs(journal, plan, expected):
    """Publish staged bytes exclusively; resume accepts only identical owned files."""
    allowed = {str(Path(path).absolute()) for path in expected}
    if {row['destination'] for row in plan} != allowed or len(plan) != len(allowed):
        raise ValueError('Checkpoint output plan does not match requested outputs')
    for row in plan:
        source = journal.path / row['source']
        if (source.is_symlink() or not source.resolve().is_relative_to(journal.path.resolve())
                or hashlib.sha256(source.read_bytes()).hexdigest() != row['sha256']):
            raise ValueError('Checkpoint staged output changed or is invalid')
        destination = Path(row['destination'])
        if destination.is_symlink():
            raise FileExistsError(f'Output is a symlink: {destination}')
        if destination.exists():
            if (not destination.is_file()
                    or hashlib.sha256(destination.read_bytes()).hexdigest() != row['sha256']):
                raise FileExistsError(f'Output differs from checkpoint: {destination}')
        else:
            # Do not hard-link the authoritative staged file to a user-editable
            # result: editing the result must not mutate checkpoint evidence.
            fd, pending = tempfile.mkstemp(prefix='.publish.', suffix='.pending', dir=journal.path)
            with os.fdopen(fd, 'wb') as stream, source.open('rb') as original:
                shutil.copyfileobj(original, stream)
                stream.flush()
                os.fsync(stream.fileno())
            os.link(pending, destination)
            _fsync_directory(destination.parent)


def output_plan(journal, staged, destinations):
    """Hash staged files and pair them with their destinations for exclusive publication."""
    rows = []
    for source, destination in zip(staged, destinations):
        with source.open('rb') as stream:
            digest = hashlib.sha256(stream.read()).hexdigest()
            os.fsync(stream.fileno())
        _fsync_directory(source.parent)
        rows.append({'source': str(source.relative_to(journal.path)),
                     'destination': str(destination.absolute()), 'sha256': digest})
    return rows
