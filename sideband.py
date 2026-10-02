"""Exact two-site NH CEST propagation using optimalcontrol (OC) operators.

All frequencies are Hz, durations seconds, relaxation/exchange rates s^-1.
N is tensor factor 0. The real orthonormal product basis has 16 components/site.
Initial state: pA*Nz, pB*Nz; observable: Nz(A)/pA. As in ONEST Matrix, no
equilibrium recovery is included (the longitudinal deviation decays to zero).
Every constant pulse segment includes RF, offsets, J, relaxation AND exchange.
"""

from functools import lru_cache

import numpy as np
from optimalcontrol.operators import E, Ix, Iy, Iz, liouvillian_comm, vec
from scipy.linalg import expm

_PAULI = (E(), 2 * Ix(), 2 * Iy(), 2 * Iz())
_BASIS = [np.kron(n, h) / 2 for n in _PAULI for h in _PAULI]
_V = np.column_stack([vec(b) for b in _BASIS])


def _generator(h):
    return (_V.conj().T @ liouvillian_comm(2 * np.pi * h) @ _V).real


_NX, _NZ = _generator(_BASIS[4]), _generator(_BASIS[12])
_HX, _HY, _HZ = (_generator(_BASIS[i]) for i in (1, 2, 3))
_J = _generator(_BASIS[15] / 2)  # Nz Hz, not 2NzHz


def composite_segments(p90_s=70e-6, cycle="RR", b1_scale=1.0):
    """Rows (duration_s, x_Hz, y_Hz); pulse lengths stay fixed as B1 changes.

    RR is the stock repeated R=90x-240y-90x; RRbar has phases +180 on Rbar.
    """
    if not np.isfinite([p90_s, b1_scale]).all() or p90_s <= 0 or b1_scale <= 0:
        raise ValueError("p90_s and b1_scale must be finite and positive")
    phases = {"R": [1], "RR": [1, 1], "RRbar": [1, -1], "MLEV4": [1, 1, -1, -1]}
    if cycle not in phases:
        raise ValueError(f"Unknown decoupling cycle: {cycle}")
    rf = b1_scale / (4 * p90_s)
    return np.array(
        [
            (d, s * x, s * y)
            for s in phases[cycle]
            for d, x, y in [(p90_s, rf, 0), (p90_s * 8 / 3, 0, rf), (p90_s, rf, 0)]
        ],
        dtype=float,
    )


def waveform_segments(path, rf_hz=None, b1_scale=1.0):
    """Read OC export_json: x/y channels, uniform time grid, metadata.pulse_dt.

    units='Hz' is absolute; units='a.u.' requires an explicit rf_hz multiplier.
    The file must contain one FULL period, including any desired supercycle.
    """
    from optimalcontrol.io import import_json

    w = import_json(path)
    if sorted(w.channels) != ["x", "y"]:
        raise ValueError("OC waveform must contain exactly x and y channels")
    dt = float(w.metadata.get("pulse_dt", np.nan))
    if (
        not np.isfinite(dt)
        or dt <= 0
        or not np.allclose(
            w.times, np.arange(len(w.times)) * dt, atol=1e-12, rtol=1e-10
        )
    ):
        raise ValueError(
            "OC waveform needs uniform times from zero and metadata.pulse_dt > 0"
        )
    if w.units == "Hz":
        if rf_hz is not None:
            raise ValueError("Do not supply rf_hz for an OC waveform already in Hz")
        scale = 1.0
    elif w.units == "a.u." and rf_hz is not None:
        scale = float(rf_hz)
    else:
        raise ValueError("OC waveform units must be Hz, or a.u. with explicit rf_hz")
    if not np.isfinite([scale, b1_scale]).all() or min(scale, b1_scale) <= 0:
        raise ValueError("Waveform RF scale must be finite and positive")
    xy = w.data[[w.channels.index("x"), w.channels.index("y")]].T * scale * b1_scale
    # Merge only exactly equal neighbours; preserve every change in the OC shape.
    starts = np.r_[0, 1 + np.flatnonzero(np.any(xy[1:] != xy[:-1], axis=1))]
    return np.column_stack((np.diff(np.r_[starts, len(xy)]) * dt, xy[starts]))


def profile(
    segments,
    offsets,
    *,
    T,
    nu,
    kab,
    kba,
    dw,
    r1,
    r2a,
    r2b,
    ha,
    hb,
    J=92.0,
    r1h=2.0,
    r2h=25.0,
):
    """I/I0 at RF offsets relative to site A; dw = nu_B - nu_A (Hz).

    ha/hb are each site's proton offset from the decoupler carrier, in Hz.
    The final fractional segment is propagated exactly, with no time rounding.
    """
    seg = np.asarray(segments, dtype=float)
    off = np.atleast_1d(np.asarray(offsets, dtype=float))
    rates = [T, nu, kab, kba, r1, r2a, r2b, r1h, r2h]
    if (
        seg.ndim != 2
        or seg.shape[1] != 3
        or len(seg) == 0
        or not np.isfinite(seg).all()
        or np.any(seg[:, 0] <= 0)
    ):
        raise ValueError("segments must be finite (n,3) rows with positive durations")
    if (
        off.ndim != 1
        or not off.size
        or not np.isfinite(off).all()
        or not np.isfinite(rates + [dw, ha, hb, J]).all()
        or min(rates) < 0
        or kba <= 0
    ):
        raise ValueError("Invalid offsets/rates; rates must be nonnegative and kba > 0")
    pa = kba / (kab + kba)
    c0 = np.zeros(32)
    c0[12], c0[28] = pa, 1 - pa
    period = seg[:, 0].sum()
    n_full = int(np.floor(T / period))
    remaining = T - n_full * period
    identity = np.eye(32)
    hrf = [np.kron(np.eye(2), a) for a in (_HX, _HY)]
    values = []
    for start in range(0, off.size, 128):
        o = off[start : start + 128]
        a = np.zeros((len(o), 32, 32))
        for j, (shift, dh, r2, rate) in enumerate(
            [(0.0, ha, r2a, kab), (dw, hb, r2b, kba)]
        ):
            blk = slice(16 * j, 16 * (j + 1))
            relaxation = (
                np.array([0, r2, r2, r1])[:, None]
                + np.array([0, r2h, r2h, r1h])[None, :]
            ).ravel()
            a[:, blk, blk] = (
                (shift - o)[:, None, None] * _NZ
                + nu * _NX
                + dh * _HZ
                + J * _J
                - np.diag(relaxation + rate)
            )
        a[:, :16, 16:], a[:, 16:, :16] = kba * np.eye(16), kab * np.eye(16)
        full = np.broadcast_to(identity, a.shape).copy()
        partial = full.copy()
        rem = remaining

        # Only identical constant segments reuse exponentials; no approximate cache.
        @lru_cache(maxsize=16)
        def propagate(duration, fx, fy):
            return expm((a + fx * hrf[0] + fy * hrf[1]) * duration)

        for duration, fx, fy in seg:
            u = propagate(duration, fx, fy)
            full = u @ full
            use = min(duration, rem)
            if use > 0:
                q = u if use == duration else propagate(use, fx, fy)
                partial = q @ partial
                rem = max(0.0, rem - use)
        end = (partial @ np.linalg.matrix_power(full, n_full)) @ c0
        values.extend(end[:, 12] / pa)
    out = np.asarray(values)
    if not np.isfinite(out).all():
        raise ValueError("Non-finite Sideband propagation")
    return out if np.ndim(offsets) else float(out[0])
