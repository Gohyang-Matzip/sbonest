"""Convert Bruker processed pseudo-2D CEST data into SBONEST text input.

Only the standard library and NumPy are used. A processed directory
(``.../pdata/1``) holds ``procs``/``proc2s`` parameter files and the ``2rr``
real part stored in submatrix blocks. Each row is one saturation offset; the
detected dimension is the proton axis. Peak intensities are read at supplied
proton positions, normalized by a reference row, and given an absolute error
from a signal-free noise region of the reference row. The saturation offsets
come from a list file (one value per row, ppm on the saturated nucleus or Hz
relative to a carrier). The converter does not guess acquisition parameters:
saturation time, nominal RF amplitude, peaks and offsets are explicit inputs.
"""
import argparse
import json
from pathlib import Path
import re

import numpy as np

_PARAM = re.compile(r"^##\$(\w+)=\s*(.*)$")
GAMMA_RATIO = {"15N": 0.101329118, "13C": 0.251449530}


def read_parameters(path):
    """Parse a Bruker JCAMP-DX parameter file into a dict of scalars and lists."""
    values = {}
    lines = Path(path).read_text(encoding="latin-1").splitlines()
    index = 0
    while index < len(lines):
        match = _PARAM.match(lines[index])
        index += 1
        if not match:
            continue
        key, raw = match.groups()
        raw = raw.strip()
        if raw.startswith("("):
            items = []
            while index < len(lines) and not lines[index].startswith("##"):
                items.extend(lines[index].split())
                index += 1
            values[key] = [_scalar(item) for item in items]
        else:
            values[key] = _scalar(raw)
    return values


def _scalar(text):
    text = text.strip()
    if text.startswith("<") and text.endswith(">"):
        return text[1:-1]
    try:
        return int(text)
    except ValueError:
        try:
            return float(text)
        except ValueError:
            return text


def read_pseudo2d(pdata):
    """Return (rows x points) real data and axis information from a pdata directory."""
    pdata = Path(pdata)
    procs = read_parameters(pdata / "procs")
    proc2s = read_parameters(pdata / "proc2s")
    for key in ("SI", "XDIM", "BYTORDP", "DTYPP", "NC_proc", "SW_p", "SF", "OFFSET"):
        if key not in procs:
            raise ValueError(f"procs lacks {key}")
    for key in ("SI", "XDIM"):
        if key not in proc2s:
            raise ValueError(f"proc2s lacks {key}")
    si2, si1 = int(procs["SI"]), int(proc2s["SI"])
    xdim2, xdim1 = int(procs["XDIM"]) or si2, int(proc2s["XDIM"]) or si1
    if si2 <= 0 or si1 <= 0 or si2 % xdim2 or si1 % xdim1:
        raise ValueError("SI must be positive multiples of XDIM in both dimensions")
    dtype = {0: np.int32, 2: np.float64}.get(int(procs["DTYPP"]))
    if dtype is None:
        raise ValueError("Unsupported DTYPP; expected 0 (int32) or 2 (float64)")
    order = ">" if int(procs["BYTORDP"]) == 1 else "<"
    raw = np.fromfile(pdata / "2rr", dtype=np.dtype(dtype).newbyteorder(order))
    if raw.size != si1 * si2:
        raise ValueError(f"2rr holds {raw.size} values; expected {si1 * si2}")
    # Submatrix storage: blocks of xdim1 rows x xdim2 points, row-major over blocks.
    blocks = raw.reshape(si1 // xdim1, si2 // xdim2, xdim1, xdim2)
    data = blocks.transpose(0, 2, 1, 3).reshape(si1, si2).astype(float)
    data *= 2.0 ** int(procs["NC_proc"])
    sf, sw, offset = float(procs["SF"]), float(procs["SW_p"]), float(procs["OFFSET"])
    ppm = offset - np.arange(si2) * (sw / sf) / si2
    return data, {"ppm": ppm, "sf_mhz": sf, "sw_hz": sw, "offset_ppm": offset, "rows": si1, "points": si2}


def read_offsets(path, *, unit, carrier_ppm=None, field_mhz=None):
    """One saturation offset per row from a text list; Hz values need a carrier.

    Bruker frequency lists (``fq1list``, ``fq2list``, ...) are accepted: a first
    line such as ``bf ppm``, ``sfo hz`` or ``P`` selects the unit of the values that
    follow, and ``O1``/``O2`` lines are ignored. An explicit ``unit`` wins only
    when the file carries no unit line.
    """
    values = []
    file_unit = None
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        text = line.split("#", 1)[0].strip()
        if not text:
            continue
        tokens = text.replace(",", " ").split()
        head = tokens[0].lower()
        if head in ("bf", "sfo", "p", "o1", "o2", "o3"):
            # Bruker list header: unit keyword, optionally followed by values on the same line.
            if head in ("bf", "sfo", "p"):
                if len(tokens) > 1 and tokens[1].lower() in ("ppm", "hz"):
                    file_unit = tokens[1].lower()
                    tokens = tokens[2:]
                else:
                    file_unit = "ppm" if head == "p" else "hz"
                    tokens = tokens[1:]
            else:
                continue
        for item in tokens:
            try:
                values.append(float(item))
            except ValueError as exc:
                raise ValueError(f"Offset list contains a non-numeric value: {item}") from exc
    if not values:
        raise ValueError("Offset list is empty")
    offsets = np.asarray(values, dtype=float)
    unit = file_unit or unit
    if unit == "ppm":
        return offsets
    if unit != "hz":
        raise ValueError("Offset unit must be ppm or hz")
    if carrier_ppm is None or field_mhz is None:
        raise ValueError("Hz offsets need --carrier-ppm and the saturated-nucleus field")
    return carrier_ppm + offsets / field_mhz


def peak_intensities(data, ppm, peaks, half_width, *, mode="max"):
    """Intensity of every peak in every row within +/- half_width ppm."""
    out = {}
    for label, position in peaks.items():
        mask = np.abs(ppm - position) <= half_width
        if not mask.any():
            raise ValueError(f"No points within {half_width} ppm of {label} at {position} ppm")
        window = data[:, mask]
        out[label] = window.max(axis=1) if mode == "max" else window.sum(axis=1)
    return out


def noise_sigma(row, ppm, region):
    """Standard deviation of one spectrum row inside a signal-free ppm region."""
    lo, hi = sorted(region)
    mask = (ppm >= lo) & (ppm <= hi)
    if mask.sum() < 8:
        raise ValueError("Noise region must contain at least eight points")
    return float(np.std(row[mask], ddof=1))


def find_peaks(row, ppm, *, threshold, min_separation_ppm=0.03, limit=50):
    """Local maxima of one spectrum row above ``threshold``, as (ppm, height) sorted by height."""
    row = np.asarray(row, dtype=float)
    candidates = np.flatnonzero((row[1:-1] > row[:-2]) & (row[1:-1] >= row[2:]) & (row[1:-1] > threshold)) + 1
    peaks = []
    for index in sorted(candidates, key=lambda i: -row[i]):
        if all(abs(ppm[index] - position) >= min_separation_ppm for position, _ in peaks):
            peaks.append((float(ppm[index]), float(row[index])))
        if len(peaks) >= limit:
            break
    return peaks


def peaks_from_reference(pdata, reference_row, noise_region, *, snr=10.0, min_separation_ppm=0.03, limit=50):
    """Candidate proton peaks (label, ppm, height) from the reference row of a pseudo-2D."""
    data, axis = read_pseudo2d(pdata)
    if not 0 <= reference_row < axis["rows"]:
        raise ValueError("Reference row outside the data")
    sigma = noise_sigma(data[reference_row], axis["ppm"], noise_region)
    found = find_peaks(data[reference_row], axis["ppm"], threshold=snr * sigma,
                       min_separation_ppm=min_separation_ppm, limit=limit)
    return [(f"P{index + 1}", position, height) for index, (position, height) in enumerate(found)]


def qa_pdf(path, data, axis, peaks, half_width, reference_row, noise_region, offsets, keep, heights, sigma):
    """QA figure: reference row with windows and noise region, then each extracted profile."""
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    ppm = axis["ppm"]
    with PdfPages(path) as pdf:
        fig, ax = plt.subplots(figsize=(11, 4.5))
        ax.plot(ppm, data[reference_row], linewidth=.7, color="black", label=f"reference row {reference_row}")
        lo, hi = sorted(noise_region)
        ax.axvspan(lo, hi, color="gray", alpha=.25, label="noise region")
        for label, position in peaks.items():
            ax.axvspan(position - half_width, position + half_width, color="tab:orange", alpha=.3)
            ax.annotate(label, (position, float(data[reference_row][np.argmin(np.abs(ppm - position))])),
                        textcoords="offset points", xytext=(0, 6), ha="center", fontsize=8)
        ax.invert_xaxis()
        ax.set(xlabel="1H chemical shift (ppm)", ylabel="intensity", title="Reference spectrum, peak windows and noise region")
        ax.legend(fontsize=8)
        pdf.savefig(fig)
        plt.close(fig)
        fig, axes = plt.subplots(1, len(peaks), figsize=(5 * len(peaks), 4), squeeze=False)
        for ax, (label, position) in zip(axes[0], peaks.items()):
            reference = heights[label][reference_row]
            ax.errorbar(offsets[keep], heights[label][keep] / reference, yerr=sigma / reference, fmt="o", markersize=3)
            ax.set(xlabel="saturation offset (ppm)", ylabel="I/I0", title=f"{label} at {position:.3f} ppm")
            ax.grid(alpha=.25)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)


def convert(pdata, offsets_path, out, *, peaks, half_width, reference_row, noise_region,
            saturation_s, v1_hz, v1err_hz=0.0, offset_unit="ppm", carrier_ppm=None,
            nucleus="15N", field_mhz=None, mode="max", dw_ppm=None, r2a=10.0, r2b=20.0,
            exclude_rows=(), qa_path=None):
    """Write one SBONEST dataset; return a summary dict."""
    data, axis = read_pseudo2d(pdata)
    if field_mhz is None:
        if nucleus not in GAMMA_RATIO:
            raise ValueError("Unknown nucleus; supply --field-mhz")
        field_mhz = axis["sf_mhz"] * GAMMA_RATIO[nucleus]
    offsets = read_offsets(offsets_path, unit=offset_unit, carrier_ppm=carrier_ppm, field_mhz=field_mhz)
    if len(offsets) != axis["rows"]:
        raise ValueError(f"Offset list has {len(offsets)} rows but the data has {axis['rows']}")
    if not 0 <= reference_row < axis["rows"]:
        raise ValueError("Reference row outside the data")
    if saturation_s <= 0 or v1_hz <= 0 or v1err_hz < 0 or half_width <= 0:
        raise ValueError("Saturation time, RF amplitude and half width must be positive; v1err nonnegative")
    heights = peak_intensities(data, axis["ppm"], peaks, half_width, mode=mode)
    sigma = noise_sigma(data[reference_row], axis["ppm"], noise_region)
    if mode == "sum":
        sigma *= np.sqrt(int(np.sum(np.abs(axis["ppm"] - next(iter(peaks.values()))) <= half_width)))
    excluded = {reference_row, *exclude_rows}
    keep = [i for i in range(axis["rows"]) if i not in excluded]
    if not keep:
        raise ValueError("No rows remain after excluding the reference")
    lines = [f"{field_mhz:.10g}", f"{saturation_s:.10g}", f"{v1_hz:.10g} {v1err_hz:.10g}",
             "# offset(ppm) intensity error  (normalized to the reference row)"]
    summary = {"pdata": str(Path(pdata).resolve()), "rows": axis["rows"], "points": axis["points"],
               "field_mhz": field_mhz, "proton_sf_mhz": axis["sf_mhz"], "reference_row": reference_row,
               "noise_sigma_raw": sigma, "peaks": {}, "excluded_rows": sorted(excluded)}
    for label, position in peaks.items():
        reference = float(heights[label][reference_row])
        if reference <= 0:
            raise ValueError(f"Reference intensity of {label} is not positive")
        normalized = heights[label][keep] / reference
        error = sigma / reference
        dw = dw_ppm.get(label, 0.0) if isinstance(dw_ppm, dict) else (dw_ppm or 0.0)
        lines.append(f"# {label} R2a: {r2a:g} R2b: {r2b:g} dw: {dw:g}")
        lines.extend(f"{o:.8g} {v:.10g} {error:.6g}" for o, v in zip(offsets[keep], normalized))
        summary["peaks"][label] = {"ppm": position, "reference_intensity": reference,
                                   "normalized_error": error, "n_points": len(keep)}
    out = Path(out)
    if out.exists() or out.is_symlink():
        raise FileExistsError(f"Output already exists: {out}")
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    summary["output"] = str(out.resolve())
    if qa_path is not None:
        qa_path = Path(qa_path)
        if qa_path.exists() or qa_path.is_symlink():
            raise FileExistsError(f"QA output already exists: {qa_path}")
        qa_pdf(qa_path, data, axis, peaks, half_width, reference_row, noise_region, offsets, keep, heights, sigma)
        summary["qa_pdf"] = str(qa_path.resolve())
    return summary


def _peak_argument(text):
    if "=" not in text:
        raise argparse.ArgumentTypeError("peaks are LABEL=ppm or LABEL=ppm:dw_ppm")
    label, rest = text.split("=", 1)
    parts = rest.split(":")
    try:
        position = float(parts[0])
        dw = float(parts[1]) if len(parts) > 1 else 0.0
    except ValueError as exc:
        raise argparse.ArgumentTypeError("peak position and dw must be numbers") from exc
    if not re.fullmatch(r"\w+\d+", label):
        raise argparse.ArgumentTypeError("peak labels must end in a residue number, e.g. A1")
    return label, position, dw


def main(argv=None):
    """Command line of the Bruker pseudo-2D converter."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("pdata", help="Bruker processed directory containing procs, proc2s and 2rr")
    parser.add_argument("offsets", help="Text list with one saturation offset per row")
    parser.add_argument("--out", required=True, help="New SBONEST dataset file")
    parser.add_argument("--peak", action="append", type=_peak_argument,
                        metavar="LABEL=ppm[:dw_ppm]", help="Proton peak position (repeatable)")
    parser.add_argument("--peaks-from-reference", action="store_true",
                        help="Without --peak: take peaks above --peak-snr times the noise from the reference row (labels P1, P2, ...)")
    parser.add_argument("--peak-snr", type=float, default=10.0, help="Signal-to-noise threshold for --peaks-from-reference")
    parser.add_argument("--qa-pdf", help="Write a QA figure (reference row, windows, noise region, profiles) to this new file")
    parser.add_argument("--half-width", type=float, default=0.02, help="Window half width in ppm")
    parser.add_argument("--mode", choices=("max", "sum"), default="max", help="Window statistic")
    parser.add_argument("--reference-row", type=int, required=True, help="Row index of the reference spectrum")
    parser.add_argument("--noise-region", type=float, nargs=2, required=True, metavar=("PPM_LO", "PPM_HI"))
    parser.add_argument("--saturation-s", type=float, required=True)
    parser.add_argument("--v1-hz", type=float, required=True)
    parser.add_argument("--v1err-hz", type=float, default=0.0)
    parser.add_argument("--offset-unit", choices=("ppm", "hz"), default="ppm")
    parser.add_argument("--carrier-ppm", type=float, help="Carrier for Hz offsets (saturated nucleus)")
    parser.add_argument("--nucleus", default="15N", choices=sorted(GAMMA_RATIO))
    parser.add_argument("--field-mhz", type=float, help="Saturated-nucleus Larmor frequency (default from SF)")
    parser.add_argument("--r2a", type=float, default=10.0)
    parser.add_argument("--r2b", type=float, default=20.0)
    parser.add_argument("--exclude-row", type=int, action="append", default=[])
    args = parser.parse_args(argv)
    if args.peak:
        peaks = {label: position for label, position, _ in args.peak}
        dw = {label: value for label, _, value in args.peak}
    elif args.peaks_from_reference:
        try:
            found = peaks_from_reference(args.pdata, args.reference_row, args.noise_region, snr=args.peak_snr)
        except (ValueError, OSError) as exc:
            parser.exit(1, f"Error: {exc}\n")
        if not found:
            parser.exit(1, "Error: no peaks above the threshold in the reference row\n")
        peaks = {label: position for label, position, _ in found}
        dw = {}
    else:
        parser.error("supply --peak LABEL=ppm or --peaks-from-reference")
    try:
        summary = convert(args.pdata, args.offsets, args.out, peaks=peaks, half_width=args.half_width,
                          reference_row=args.reference_row, noise_region=args.noise_region,
                          saturation_s=args.saturation_s, v1_hz=args.v1_hz, v1err_hz=args.v1err_hz,
                          offset_unit=args.offset_unit, carrier_ppm=args.carrier_ppm, nucleus=args.nucleus,
                          field_mhz=args.field_mhz, mode=args.mode, dw_ppm=dw, r2a=args.r2a, r2b=args.r2b,
                          exclude_rows=args.exclude_row, qa_path=args.qa_pdf)
    except (ValueError, OSError) as exc:
        parser.exit(1, f"Error: {exc}\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
