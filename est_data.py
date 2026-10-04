"""
Original 2015 by Donghan Lee
Python 3 conversion: 2025
"""

import re
import numpy as np


_NUMBER = r"([+-]?\d+\.?\d*(?:[eE][+-]?\d+)?)"
_SCALAR_RE = re.compile(rf"^\s*{_NUMBER}\s*(?:#.*)?$")
_RF_RE = re.compile(rf"^\s*{_NUMBER}\s+{_NUMBER}\s*(?:#.*)?$")
_DATA_RE = re.compile(rf"^\s*{_NUMBER}\s+{_NUMBER}\s+{_NUMBER}\s*(?:#.*)?$")
_RESIDUE_SIMPLE_RE = re.compile(r"^#\s*(\w+\d+)\s*(?:#.*)?$")
_RESIDUE_FULL_RE = re.compile(
    rf"^#\s*(\w+\d+)\s+R2a:\s*{_NUMBER}\s+R2b:\s*{_NUMBER}\s+dw:\s*{_NUMBER}\s*(?:#.*)?$"
)


def _read_conditions(stream, file_name):
    """Read field, saturation time, RF amplitude, and RF uncertainty in order."""
    values = []
    for pattern, label, error_label in (
        (_SCALAR_RE, "B0 field", "Field"),
        (_SCALAR_RE, "T", "T"),
        (_RF_RE, "V1", "V1"),
    ):
        line = stream.readline()
        if not line:
            raise ValueError(f"Unexpected EOF at {label} line in {file_name}")
        match = pattern.match(line)
        if not match:
            raise ValueError(f"{error_label} line parse error: '{line.strip()}'")
        values.extend(map(float, match.groups()))
    return values


def _read_points(stream, spectrum, label, file_name, add_error):
    """Fill one spectrum and return the next comment/header, or an empty EOF."""
    for line in stream:
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("#"):
            break
        match = _DATA_RE.match(line)
        if not match:
            raise ValueError(
                f"Malformed data line in {file_name} under residue {label}: {stripped}"
            )
        offset, intensity, std = map(float, match.groups())
        if add_error:
            intensity += std * np.random.randn()
        spectrum.offset.append(offset)
        spectrum.int.append(intensity)
        spectrum.intstd.append(std)
    else:
        line = ""
    if not spectrum.offset:
        raise ValueError(f"No data points for residue {label}")
    return line


class EstSpec:
    def __init__(self):
        self.field = 0.0
        self.T = 0.0
        self.centerppm = 0.0
        self.v1 = 0.0
        self.v1err = 0.0
        self.offset = []
        self.int = []
        self.intstd = []
        self.initdw = 0.1
        self.initr2a = 10.0
        self.initr2b = 100.0

    def info(self):
        print(f"--- B0 field: {self.field:8.3f} [MHz]")
        print(f"-          T: {1000.0 * self.T:8.3f} [ms]")
        print(f"-         v1: {self.v1:8.3f} {self.v1err:8.3f} [Hz]")
        for offset_val, intensity_val, std_val in zip(
            self.offset, self.int, self.intstd
        ):
            print(f"{offset_val:8.3f}  {intensity_val:8.3f}   {std_val:8.3f}")
        print("---")


class Residue:
    def __init__(self):
        self.label = ""
        self.estSpecs = []
        self.active = True


class EstDataSet:
    def __init__(self):
        self.res = []
        self.fields = []
        self.Ts = []
        self.centerppms = []
        self.v1s = []
        self.v1errs = []
        self.initR2 = False

    def addData(self, fileName, add_error_to_intensity=False, add_error_to_v1=False):
        """Append ONEST spectra, retaining residue and per-file acquisition order."""
        try:
            with open(fileName, "r") as stream:
                field, duration, v1, v1err = _read_conditions(stream, fileName)
                if add_error_to_v1:
                    v1 += v1err * np.random.randn()
                self.fields.append(field)
                self.Ts.append(duration)
                self.v1s.append(v1)
                self.v1errs.append(v1err)

                header = stream.readline()
                if not header:
                    raise ValueError(
                        f"Unexpected EOF after V1 line (expected header) in {fileName}"
                    )
                if (
                    _RESIDUE_SIMPLE_RE.match(header.strip())
                    or _RESIDUE_FULL_RE.match(header.strip())
                    or _DATA_RE.match(header)
                ):
                    raise ValueError(
                        "Missing column-header line after V1; found a residue or data row"
                    )

                file_has_initial = None
                line = stream.readline()
                while line:
                    stripped = line.strip()
                    match = _RESIDUE_FULL_RE.match(line)
                    has_initial = match is not None
                    if match is None:
                        match = _RESIDUE_SIMPLE_RE.match(line)
                    if not stripped or (stripped.startswith("#") and match is None):
                        line = stream.readline()
                        continue
                    if match is None:
                        raise ValueError(
                            "Data or unexpected text outside a residue block; "
                            "check residue headers and comments inside data blocks: "
                            f"{stripped}"
                        )
                    if file_has_initial is None:
                        file_has_initial = has_initial
                    elif file_has_initial != has_initial:
                        expected = "expected R2a/b" if has_initial else "no R2a/b expected"
                        raise ValueError(
                            f"Inconsistent residue format in {fileName} ({expected}). Line: {stripped}"
                        )

                    label = match.group(1)
                    residue = next((r for r in self.res if r.label == label), None)
                    if residue is None:
                        residue = Residue()
                        residue.label = label
                        self.res.append(residue)
                    spectrum = EstSpec()
                    spectrum.field, spectrum.T, spectrum.v1, spectrum.v1err = (
                        field, duration, v1, v1err
                    )
                    if has_initial:
                        spectrum.initr2a, spectrum.initr2b, spectrum.initdw = map(
                            float, match.groups()[1:]
                        )
                    line = _read_points(
                        stream, spectrum, label, fileName, add_error_to_intensity
                    )
                    residue.estSpecs.append(spectrum)

                if file_has_initial is None:
                    raise ValueError("No residue data found")
                if self.initR2 is False and file_has_initial is True:
                    self.initR2 = True
        except ValueError as e:
            raise ValueError(f"Error processing file {fileName}: {e}") from e

    def addDataWithError(self, fileName):
        self.addData(fileName, add_error_to_intensity=True)

    def getResidues(self):
        return [{"name": r.label, "flag": "on"} for r in self.res]

    def info(self):
        for i, (field, T_val, v1_val, v1err_val) in enumerate(
            zip(self.fields, self.Ts, self.v1s, self.v1errs)
        ):
            print(f"\nExperimental Condition Set {i + 1}:")
            print(f"  B0 Field: {field:8.3f} [MHz]")
            print(f"  T:        {1000.0 * T_val:8.3f} [ms]")
            print(f"  v1:       {v1_val:8.3f} +/- {v1err_val:8.3f} [Hz]")
            print(
                f"{'#Residue':<10} {'Offset (ppm)':<15} {'Intensity':<12} {'StdDev':<10}"
            )
            print("-" * 50)

            for r_obj in self.res:
                for ep in r_obj.estSpecs:
                    if (
                        ep.field == field
                        and ep.T == T_val
                        and ep.v1 == v1_val
                        and ep.v1err == v1err_val
                    ):
                        for j_idx in range(len(ep.offset)):
                            print(
                                f"{r_obj.label:<10} {ep.offset[j_idx]:<15.3f} {ep.int[j_idx]:<12.3f} {ep.intstd[j_idx]:<10.3f}"
                            )
