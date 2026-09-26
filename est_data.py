"""
Original 2015 by Donghan Lee
Python 3 conversion: 2025
"""

import re
import numpy as np


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

    def _parse_and_check(
        self, line, regex_pattern, expected_groups, error_message_prefix
    ):
        # Attempt to strip comments first for lines that might have them universally
        # For V1_RE, the comment handling is now part of the regex itself.
        # For other lines, if they can also have comments, this would need adjustment.
        line_to_parse = (
            line  # For FIELD_RE and T_RE, comments are not expected by original parser.
        )

        match = regex_pattern.match(line_to_parse)
        if not match or len(match.groups()) != expected_groups:
            # For debugging, show what the regex captured vs expected
            # captured_groups = match.groups() if match else "No match"
            # print(f"Debug: Line='{line.strip()}', Regex='{regex_pattern.pattern}', ExpectedGroups={expected_groups}, Captured={captured_groups}")
            raise ValueError(f"{error_message_prefix}: '{line.strip()}'")
        return match.groups()

    def addData(self, fileName, add_error_to_intensity=False, add_error_to_v1=False):
        FIELD_RE = re.compile(
            r"^\s*([+-]?\d+\.?\d*(?:[eE][+-]?\d+)?)\s*(?:#.*)?$"
        )  # Allow comment
        T_RE = re.compile(
            r"^\s*([+-]?\d+\.?\d*(?:[eE][+-]?\d+)?)\s*(?:#.*)?$"
        )  # Allow comment
        # Regex for V1 line, allowing for an optional comment at the end
        V1_RE = re.compile(
            r"^\s*([+-]?\d+\.?\d*(?:[eE][+-]?\d+)?)\s+([+-]?\d+\.?\d*(?:[eE][+-]?\d+)?)\s*(?:#.*)?$"
        )
        RESIDUE_SIMPLE_RE = re.compile(r"^#\s*(\w+\d+)\s*(?:#.*)?$")  # Allow comment
        RESIDUE_FULL_RE = re.compile(
            r"^#\s*(\w+\d+)\s+R2a:\s*([+-]?\d+\.?\d*(?:[eE][+-]?\d+)?)\s+R2b:\s*([+-]?\d+\.?\d*(?:[eE][+-]?\d+)?)\s+dw:\s*([+-]?\d+\.?\d*(?:[eE][+-]?\d+)?)\s*(?:#.*)?$"
        )  # Allow comment
        DATA_LINE_RE = re.compile(
            r"^\s*([+-]?\d+\.?\d*(?:[eE][+-]?\d+)?)\s+([+-]?\d+\.?\d*(?:[eE][+-]?\d+)?)\s+([+-]?\d+\.?\d*(?:[eE][+-]?\d+)?)\s*(?:#.*)?$"
        )  # Allow comment for data lines too

        try:
            with open(fileName, "r") as inFile:
                # Read B0 field
                line1 = inFile.readline()
                if not line1:
                    raise ValueError(f"Unexpected EOF at B0 field line in {fileName}")
                currentField = float(
                    self._parse_and_check(line1, FIELD_RE, 1, "Field line parse error")[
                        0
                    ]
                )

                # Read T
                line2 = inFile.readline()
                if not line2:
                    raise ValueError(f"Unexpected EOF at T line in {fileName}")
                currentT_val = float(
                    self._parse_and_check(line2, T_RE, 1, "T line parse error")[0]
                )

                # Read v1
                line3 = inFile.readline()
                if not line3:
                    raise ValueError(f"Unexpected EOF at V1 line in {fileName}")
                v1_groups = self._parse_and_check(
                    line3, V1_RE, 2, "V1 line parse error"
                )
                currentV1_val = float(v1_groups[0])
                currentV1err_val = float(v1_groups[1])

                if add_error_to_v1:
                    currentV1_val += currentV1err_val * np.random.randn()

                self.fields.append(currentField)
                self.Ts.append(currentT_val)
                self.v1s.append(currentV1_val)
                self.v1errs.append(currentV1err_val)

                # Skip the fourth line (usually a header/comment like "#offset(ppm) Intensity error")
                line4 = inFile.readline()
                if not line4:
                    raise ValueError(
                        f"Unexpected EOF after V1 line (expected header) in {fileName}"
                    )

                current_file_has_initR2_info = None

                line = inFile.readline()
                while line:
                    line_stripped = (
                        line.strip()
                    )  # Use stripped line for checks, original line for parsing
                    if not line_stripped:  # Skip blank lines
                        line = inFile.readline()
                        continue

                    resLabel, initR2a, initR2b, initDW = None, None, None, None
                    is_residue_line = False

                    if line_stripped.startswith(
                        "#"
                    ):  # All residue lines must start with #
                        match_full = RESIDUE_FULL_RE.match(
                            line
                        )  # Try full format first
                        if match_full:
                            if current_file_has_initR2_info is None:
                                current_file_has_initR2_info = True
                            if not current_file_has_initR2_info:
                                raise ValueError(
                                    f"Inconsistent residue format in {fileName} (expected R2a/b). Line: {line_stripped}"
                                )
                            resLabel, r2a_str, r2b_str, dw_str = match_full.groups()[
                                :4
                            ]  # Regex has 4 groups before optional comment
                            initR2a, initR2b, initDW = (
                                float(r2a_str),
                                float(r2b_str),
                                float(dw_str),
                            )
                            is_residue_line = True
                        else:
                            match_simple = RESIDUE_SIMPLE_RE.match(
                                line
                            )  # Try simple format
                            if match_simple:
                                if current_file_has_initR2_info is None:
                                    current_file_has_initR2_info = False
                                if current_file_has_initR2_info:
                                    raise ValueError(
                                        f"Inconsistent residue format in {fileName} (no R2a/b expected). Line: {line_stripped}"
                                    )
                                resLabel = match_simple.group(1)
                                is_residue_line = True
                            else:  # Not a recognized residue line, could be general comment. Skip.
                                line = inFile.readline()
                                continue

                    if (
                        not is_residue_line
                    ):  # Should be data or error if not residue line
                        # This path might be taken if a non-comment, non-data line appears where data is expected.
                        # Based on syn10.txt, data lines do not start with #.
                        # If it's not a residue line here, it means we are past a residue block or file is malformed.
                        # For now, if it's not a data line, it's an issue.
                        data_match_check = DATA_LINE_RE.match(line)
                        if not data_match_check:
                            # print(f"Warning in {fileName}: Unexpected line, skipping: {line_stripped}")
                            line = inFile.readline()
                            continue
                        # If it was data, it would have been consumed by the inner loop of the previous residue.
                        # This implies an issue with file structure or previous loop termination.
                        # Safest to assume an error or try to find next residue block.
                        # For now, skip. This can happen if there are blank lines between data blocks.
                        line = inFile.readline()
                        continue

                    resid_obj = next((r for r in self.res if r.label == resLabel), None)
                    if resid_obj is None:
                        resid_obj = Residue()
                        resid_obj.label = resLabel
                        self.res.append(resid_obj)

                    ep = EstSpec()
                    ep.field, ep.T, ep.v1, ep.v1err = (
                        currentField,
                        currentT_val,
                        currentV1_val,
                        currentV1err_val,
                    )
                    if current_file_has_initR2_info:
                        ep.initr2a, ep.initr2b, ep.initdw = initR2a, initR2b, initDW

                    # Read data points for this residue
                    line = (
                        inFile.readline()
                    )  # Read next line, which should be first data line or new residue
                    while line:
                        line_stripped_data = line.strip()
                        if not line_stripped_data:  # Skip blank lines within data block
                            line = inFile.readline()
                            continue
                        if line_stripped_data.startswith(
                            "#"
                        ):  # End of data or new residue
                            break

                        data_match = DATA_LINE_RE.match(line)
                        if data_match:
                            offset, intensity_val, std_val = map(
                                float, data_match.groups()[:3]
                            )  # Regex has 3 groups before comment
                            if add_error_to_intensity:
                                intensity_val += std_val * np.random.randn()
                            ep.offset.append(offset)
                            ep.int.append(intensity_val)
                            ep.intstd.append(std_val)
                        else:
                            raise ValueError(
                                f"Malformed data line in {fileName} under residue {resLabel}: {line_stripped_data}"
                            )
                        line = (
                            inFile.readline()
                        )  # Read next data line or next residue/EOF

                    if ep.offset:
                        resid_obj.estSpecs.append(ep)
                    # 'line' now holds the next residue header, a comment, a blank line, or is None (EOF)

                if self.initR2 is False and current_file_has_initR2_info is True:
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
