###
# Original 2015 by Donghan Lee
# Python 3 / scipy.optimize.least_squares port: 2025
###

import numpy as np
from scipy.linalg import expm
import concurrent.futures

from matplotlib.pyplot import figure, close
from matplotlib.backends.backend_pdf import PdfPages
from os import getlogin
from platform import uname
from time import ctime

from est_data import EstDataSet
import fit as fit_module

VALID_METHODS = [
    "Baldwin",
    "Matrix",
    "NoEx",
    "Matrix_3st_Linear",
    "Matrix_3st_Triangle",
]

# Per-method flat-parameter layout: (global keys, per-residue keys, per-residue defaults)
PARAM_LAYOUT = {
    "NoEx": ([], ["dGs", "r1s", "r2as"], [0.0, 1.0, 10.0]),
    "Baldwin": (
        ["kab", "kba"],
        ["dGs", "dws", "r1s", "r2as", "r2bs"],
        [0.0, 0.1, 1.0, 10.0, 20.0],
    ),
    "Matrix_3st_Linear": (
        ["kab", "kba", "kbc", "kcb"],
        ["dGs", "dws", "dwCs", "r1s", "r2as", "r2bs", "r2cs"],
        [0.0, 0.1, 0.2, 1.0, 10.0, 20.0, 20.0],
    ),
    "Matrix_3st_Triangle": (
        ["kab", "kba", "kbc", "kcb", "kca", "kac"],
        ["dGs", "dws", "dwCs", "r1s", "r2as", "r2bs", "r2cs"],
        [0.0, 0.1, 0.2, 1.0, 10.0, 20.0, 20.0],
    ),
}
PARAM_LAYOUT["Matrix"] = PARAM_LAYOUT["Baldwin"]

# Per-residue labels for text logs (padding matches legacy output format)
RES_LOG_LABELS = {
    "NoEx": ["peak pos [ppm]: ", "R1 [s-1]:       ", "R2a [s-1]:      "],
    "Baldwin": [
        "peak pos [ppm]: ",
        "cs diff [ppm]:  ",
        "R1 [s-1]:       ",
        "R2a [s-1]:      ",
        "R2b [s-1]:      ",
    ],
    "Matrix_3st_Linear": [
        "peak pos [ppm]: ",
        "cs diff B [ppm]:",
        "cs diff C [ppm]:",
        "R1 [s-1]:       ",
        "R2a [s-1]:      ",
        "R2b [s-1]:      ",
        "R2c [s-1]:      ",
    ],
}
RES_LOG_LABELS["Matrix"] = RES_LOG_LABELS["Baldwin"]
RES_LOG_LABELS["Matrix_3st_Triangle"] = RES_LOG_LABELS["Matrix_3st_Linear"]


def fast_gaussian(x, mu, sigma):
    """Gaussian PDF, faster than scipy.stats.norm.pdf."""
    if sigma == 0:
        return np.where(x == mu, 1.0, 0.0)
    return (1.0 / (np.sqrt(2 * np.pi) * sigma)) * np.exp(-0.5 * ((x - mu) / sigma) ** 2)


def b1_weights(v1, v1err):
    """B1 inhomogeneity sampling: (angular frequencies, normalized weights)."""
    w1 = v1 * 2.0 * np.pi
    w1err = v1err * 2.0 * np.pi
    if w1err == 0:
        return np.array([w1]), np.array([1.0])
    wxs = np.linspace(-2.0 * w1err + w1, 2.0 * w1err + w1, 10)
    weights = fast_gaussian(wxs, w1, w1err)
    s = np.sum(weights)
    if s == 0 or np.isnan(s):
        return wxs, np.ones_like(wxs) / len(wxs)
    return wxs, weights / s


def matrix_calc_worker_chunk(args):
    """
    2-state Bloch-McConnell numerical propagation over a chunk of offsets.
    Module-level so ProcessPoolExecutor can pickle it.
    args: (kab, kba, dG, dw, R1, R2a, R2b, offsets, T, v1, v1err, B0)
    """
    kab, kba, dG, dw, R1, R2a, R2b, offsets, T, v1, v1err, B0 = args

    R1a = R1b = R1
    wG = dG * B0 * 2.0 * np.pi
    wE = (dG + dw) * B0 * 2.0 * np.pi
    wy = 0.0

    kex = kab + kba
    pB, pA = (kab / kex, kba / kex) if kex != 0 else (0.0, 1.0)

    wxs, weights = b1_weights(v1, v1err)
    refM = pA if pA != 0 else 1.0
    startM = np.array([0, 0, 0, pA, 0, 0, pB], dtype=float)

    results = []
    for dRF in offsets:
        wRF = dRF * B0 * 2.0 * np.pi
        wa = wG - wRF
        wb = wE - wRF
        preIcal = 0.0
        for wx, w in zip(wxs, weights):
            # Basis: [1, MxA, MyA, MzA, MxB, MyB, MzB]
            Alist = [
                [0, 0, 0, 0, 0, 0, 0],
                [0, -kab - R2a, -wa, wy, kba, 0, 0],
                [0, wa, -kab - R2a, -wx, 0, kba, 0],
                [2 * R1a * pA, -wy, wx, -kab - R1a, 0, 0, kba],
                [0, kab, 0, 0, -kba - R2b, -wb, wy],
                [0, 0, kab, 0, wb, -kba - R2b, -wx],
                [2 * R1b * pB, 0, 0, kab, -wy, wx, -kba - R1b],
            ]
            endM = np.dot(expm(np.array(Alist, dtype=float) * T), startM)
            preIcal += w * endM[3] / refM
        results.append(max(0.0, preIcal))
    return results


class est_model:
    def __init__(self):
        self.programName = "est ver. 2.1 (Python 3, Scipy LS, Refactored Fit)"
        self.method = "Baldwin"
        self.dataset = EstDataSet()
        self.verbose = False

        self.chi2 = 0.0
        self.npar = 0
        self.dof = 0
        self.nvar = 0

        self.executor = None  # ProcessPoolExecutor for the Matrix method

    # --- Calculation models ---

    def matrix_calc_noex(self, dG, R1, R2a, dRF, T, v1, v1err, B0):
        """Vectorized analytical R1rho decay, no exchange. dRF scalar or array."""
        dRF_arr = np.atleast_1d(dRF)
        wa = (dG - dRF_arr) * B0 * 2.0 * np.pi

        wxs, weights = b1_weights(v1, v1err)
        w1_loop = wxs[:, np.newaxis]
        wa_grid = wa[np.newaxis, :]

        weff2 = wa_grid**2 + w1_loop**2
        sintheta2 = w1_loop**2 / weff2
        costheta2 = wa_grid**2 / weff2
        R1rho = R1 * costheta2 + R2a * sintheta2

        exp_arg = -1.0 * T * R1rho
        exp_val = np.exp(exp_arg)
        exp_val[exp_arg < -700] = 0.0

        Icalc = np.sum(weights[:, np.newaxis] * costheta2 * exp_val, axis=0)
        if np.ndim(dRF) == 0:
            return max(0.0, Icalc.item())
        return np.maximum(0.0, Icalc)

    def matrix_calc_3st_linear(
        self,
        kab,
        kba,
        kbc,
        kcb,
        dG,
        dw,
        dwC,
        R1,
        R2a,
        R2b,
        R2c,
        dRF,
        T,
        v1,
        v1err,
        B0,
        observable="A",
    ):
        """3-state linear exchange A <-> B <-> C. observable: 'A' or 'B'."""
        R1a = R1b = R1c = R1
        wRF = dRF * B0 * 2.0 * np.pi
        wa = dG * B0 * 2.0 * np.pi - wRF
        wb = (dG + dw) * B0 * 2.0 * np.pi - wRF
        wc = (dG + dwC) * B0 * 2.0 * np.pi - wRF
        wy = 0.0

        # Equilibrium populations from detailed balance
        ratio_ba = kab / kba if kba != 0 else 0
        ratio_cb = kbc / kcb if kcb != 0 else 0
        denom = 1.0 + ratio_ba + ratio_ba * ratio_cb
        if denom == 0:
            pA, pB, pC = 1.0, 0.0, 0.0
        else:
            pA = 1.0 / denom
            pB = ratio_ba * pA
            pC = ratio_ba * ratio_cb * pA

        wxs, weights = b1_weights(v1, v1err)
        startM = np.array([0, 0, 0, pA, 0, 0, pB, 0, 0, pC], dtype=float)
        obs_idx = 6 if observable == "B" else 3
        refM = {3: pA, 6: pB}[obs_idx] or 1.0

        preIcal = 0.0
        for wx, w in zip(wxs, weights):
            # Basis: [1, MxA, MyA, MzA, MxB, MyB, MzB, MxC, MyC, MzC]
            Alist = [
                [0] * 10,
                [0, -R2a - kab, -wa, wy, kba, 0, 0, 0, 0, 0],
                [0, wa, -R2a - kab, -wx, 0, kba, 0, 0, 0, 0],
                [2 * R1a * pA, -wy, wx, -R1a - kab, 0, 0, kba, 0, 0, 0],
                [0, kab, 0, 0, -R2b - kba - kbc, -wb, wy, kcb, 0, 0],
                [0, 0, kab, 0, wb, -R2b - kba - kbc, -wx, 0, kcb, 0],
                [2 * R1b * pB, 0, 0, kab, -wy, wx, -R1b - kba - kbc, 0, 0, kcb],
                [0, 0, 0, 0, kbc, 0, 0, -R2c - kcb, -wc, wy],
                [0, 0, 0, 0, 0, kbc, 0, wc, -R2c - kcb, -wx],
                [2 * R1c * pC, 0, 0, 0, 0, 0, kbc, -wy, wx, -R1c - kcb],
            ]
            endM = np.dot(expm(np.array(Alist, dtype=float) * T), startM)
            preIcal += w * endM[obs_idx] / refM
        return max(0.0, preIcal)

    def matrix_calc_3st_triangle(
        self,
        kab,
        kba,
        kbc,
        kcb,
        kca,
        kac,
        dG,
        dw,
        dwC,
        R1,
        R2a,
        R2b,
        R2c,
        dRF,
        T,
        v1,
        v1err,
        B0,
    ):
        """3-state triangular exchange A <-> B <-> C <-> A. A is observable."""
        R1a = R1b = R1c = R1
        wRF = dRF * B0 * 2.0 * np.pi
        wa = dG * B0 * 2.0 * np.pi - wRF
        wb = (dG + dw) * B0 * 2.0 * np.pi - wRF
        wc = (dG + dwC) * B0 * 2.0 * np.pi - wRF
        wy = 0.0

        # Equilibrium populations: solve K p = 0 with sum(p) = 1
        K_mat = np.array(
            [
                [-kab - kac, kba, kca],
                [kab, -kba - kbc, kcb],
                [1, 1, 1],
            ]
        )
        try:
            pA, pB, pC = np.linalg.solve(K_mat, np.array([0, 0, 1]))
        except np.linalg.LinAlgError:
            pA, pB, pC = 1.0, 0.0, 0.0

        wxs, weights = b1_weights(v1, v1err)
        startM = np.array([0, 0, 0, pA, 0, 0, pB, 0, 0, pC], dtype=float)
        refM = pA if pA != 0 else 1.0

        preIcal = 0.0
        for wx, w in zip(wxs, weights):
            # Basis: [1, MxA, MyA, MzA, MxB, MyB, MzB, MxC, MyC, MzC]
            Alist = [
                [0] * 10,
                [0, -R2a - kab - kac, -wa, wy, kba, 0, 0, kca, 0, 0],
                [0, wa, -R2a - kab - kac, -wx, 0, kba, 0, 0, kca, 0],
                [2 * R1a * pA, -wy, wx, -R1a - kab - kac, 0, 0, kba, 0, 0, kca],
                [0, kab, 0, 0, -R2b - kba - kbc, -wb, wy, kcb, 0, 0],
                [0, 0, kab, 0, wb, -R2b - kba - kbc, -wx, 0, kcb, 0],
                [2 * R1b * pB, 0, 0, kab, -wy, wx, -R1b - kba - kbc, 0, 0, kcb],
                [0, kac, 0, 0, kbc, 0, 0, -R2c - kcb - kca, -wc, wy],
                [0, 0, kac, 0, 0, kbc, 0, wc, -R2c - kcb - kca, -wx],
                [2 * R1c * pC, 0, 0, kac, 0, 0, kbc, -wy, wx, -R1c - kcb - kca],
            ]
            endM = np.dot(expm(np.array(Alist, dtype=float) * T), startM)
            preIcal += w * endM[3] / refM
        return max(0.0, preIcal)

    def matrix_calc(self, kab, kba, dG, dw, R1, R2a, R2b, dRF, T, v1, v1err, B0):
        """2-state numerical propagation for a single offset."""
        return matrix_calc_worker_chunk(
            (kab, kba, dG, dw, R1, R2a, R2b, [dRF], T, v1, v1err, B0)
        )[0]

    def Baldwin(self, kab, kba, dG, dw, R1, R2a, R2b, dRF, T, v1, v1err, B0):
        """Vectorized Baldwin analytical 2-state model. dRF scalar or array."""
        epsilon = 1e-9
        kex = kab + kba
        pE, pG = (kab / kex, kba / kex) if kex != 0 else (0.0, 1.0)
        dR = R2b - R2a

        dRF_arr = np.atleast_1d(dRF)
        wRF = dRF_arr * B0 * 2.0 * np.pi
        wG = dG * B0 * 2.0 * np.pi
        ddG = wG - wRF
        wE = (dG + dw) * B0 * 2.0 * np.pi
        ddE = wE - wRF
        ddw = ddE - ddG  # = dw*B0*2pi, independent of offset
        da = pG * ddG + pE * ddE

        wxs, weights = b1_weights(v1, v1err)

        # Broadcast: w1 samples (10,1) x offsets (1,N)
        w1_loop = wxs[:, np.newaxis]
        ddG_grid = ddG[np.newaxis, :]
        ddE_grid = ddE[np.newaxis, :]
        da_grid = da[np.newaxis, :]

        OG2 = w1_loop**2 + ddG_grid**2
        OE2 = w1_loop**2 + ddE_grid**2
        Oa2 = w1_loop**2 + da_grid**2

        tantheta2 = (w1_loop**2) / (da_grid**2 + epsilon)
        sintheta2 = (w1_loop**2) / (Oa2 + epsilon)
        costheta2 = (da_grid**2) / (Oa2 + epsilon)

        F1p = pG * pE * (ddw**2)
        F2p = (
            (kex**2) + (w1_loop**2) + ((ddG_grid * ddE_grid / (da_grid + epsilon)) ** 2)
        )
        Dp = (kex**2) + OG2 * OE2 / (Oa2 + epsilon)

        F1 = pE * (OG2 + (kex**2) + dR * pG * kex)
        F2 = 2 * kex + (w1_loop**2) / (kex + epsilon) + dR * pG
        F3 = 3 * pE * kex + (
            2 * pG * kex
            + (w1_loop**2) / (kex + epsilon)
            + dR
            + dR * ((pE * kex) ** 2) / (OG2 + epsilon)
        ) * (OG2 / (w1_loop**2 + epsilon))

        den = Dp + dR * F3 * sintheta2 + epsilon
        CR1 = (F2p + (F1p + dR * (F3 - F2)) * tantheta2) / den
        CR2 = (
            Dp / (sintheta2 + epsilon) - F2p / (tantheta2 + epsilon) - F1p + dR * F2
        ) / den
        Rex = (F1p * kex + dR * F1) / den

        R1r = CR1 * R1 * costheta2 + (CR2 * R2a + Rex) * sintheta2

        exp_arg = -1.0 * T * R1r
        exp_val = np.exp(exp_arg)
        exp_val[exp_arg < -700] = 0.0

        Icalc = np.sum(weights[:, np.newaxis] * costheta2 * exp_val, axis=0)
        if np.ndim(dRF) == 0:
            return max(0.0, Icalc.item())
        return np.maximum(0.0, Icalc)

    # --- Parameter handling ---

    def selMethod(self, initConf):
        """Sets the calculation method based on initConf."""
        method_val = (
            initConf.get("Method", "Baldwin") if initConf is not None else self.method
        )
        self.method = method_val if method_val in VALID_METHODS else "Baldwin"
        if self.verbose:
            print(f"Selected method: {self.method}")

    def seParam(self, p_flat):
        """Unpack flat parameter vector into {name: value or per-residue list} for self.method."""
        gkeys, rkeys, rdefaults = PARAM_LAYOUT[self.method]
        if len(p_flat) < len(gkeys):
            raise ValueError(
                f"Param array p_flat (len {len(p_flat)}) too short for {self.method}."
            )
        P = {k: p_flat[i] for i, k in enumerate(gkeys)}
        for k in rkeys:
            P[k] = []
        idx = len(gkeys)
        for res in self.dataset.res:
            if not res.active:
                for k in rkeys:
                    P[k].append(0.0)
            elif idx + len(rkeys) <= len(p_flat):
                for j, k in enumerate(rkeys):
                    P[k].append(p_flat[idx + j])
                idx += len(rkeys)
            else:
                for k, d in zip(rkeys, rdefaults):
                    P[k].append(d)
        return P

    def calc(self, P, i, dRF, es):
        """Calculated intensity for residue index i at offset(s) dRF of spectrum es. Scalar in, scalar out."""
        a = (es.T, es.v1, es.v1err, es.field)
        if self.method == "Baldwin":
            return self.Baldwin(
                P["kab"],
                P["kba"],
                P["dGs"][i],
                P["dws"][i],
                P["r1s"][i],
                P["r2as"][i],
                P["r2bs"][i],
                dRF,
                *a,
            )
        if self.method == "NoEx":
            return self.matrix_calc_noex(
                P["dGs"][i], P["r1s"][i], P["r2as"][i], dRF, *a
            )
        if self.method == "Matrix":
            vals = [
                self.matrix_calc(
                    P["kab"],
                    P["kba"],
                    P["dGs"][i],
                    P["dws"][i],
                    P["r1s"][i],
                    P["r2as"][i],
                    P["r2bs"][i],
                    o,
                    *a,
                )
                for o in np.atleast_1d(dRF)
            ]
        elif self.method == "Matrix_3st_Linear":
            vals = [
                self.matrix_calc_3st_linear(
                    P["kab"],
                    P["kba"],
                    P["kbc"],
                    P["kcb"],
                    P["dGs"][i],
                    P["dws"][i],
                    P["dwCs"][i],
                    P["r1s"][i],
                    P["r2as"][i],
                    P["r2bs"][i],
                    P["r2cs"][i],
                    o,
                    *a,
                    observable="A",
                )
                for o in np.atleast_1d(dRF)
            ]
        else:  # Matrix_3st_Triangle
            vals = [
                self.matrix_calc_3st_triangle(
                    P["kab"],
                    P["kba"],
                    P["kbc"],
                    P["kcb"],
                    P["kca"],
                    P["kac"],
                    P["dGs"][i],
                    P["dws"][i],
                    P["dwCs"][i],
                    P["r1s"][i],
                    P["r2as"][i],
                    P["r2bs"][i],
                    P["r2cs"][i],
                    o,
                    *a,
                )
                for o in np.atleast_1d(dRF)
            ]
        return np.array(vals) if np.ndim(dRF) else vals[0]

    # --- Fitting ---

    def errFunc(self, p_flat):
        try:
            P = self.seParam(p_flat)
        except ValueError as e:
            if self.verbose:
                print(f"Error in seParam (errFunc): {e}")
            n_residuals = (
                self.nvar
                if self.nvar > 0
                else max(
                    1,
                    sum(
                        len(es.offset)
                        for r in self.dataset.res
                        if r.active
                        for es in r.estSpecs
                    ),
                )
            )
            return np.full(n_residuals, 1e6, dtype=float)

        active = [(i, res) for i, res in enumerate(self.dataset.res) if res.active]

        # Chunked multiprocessing path for the (slow) 2-state Matrix method
        CHUNK_SIZE = 50
        total_offsets = sum(
            len(res.estSpecs[0].offset) for _, res in active if res.estSpecs
        )
        use_parallel = (
            self.method == "Matrix"
            and self.executor is not None
            and total_offsets > 200
        )

        if use_parallel:
            tasks = []
            for i, res in active:
                for es in res.estSpecs:
                    for k in range(0, len(es.offset), CHUNK_SIZE):
                        tasks.append(
                            (
                                P["kab"],
                                P["kba"],
                                P["dGs"][i],
                                P["dws"][i],
                                P["r1s"][i],
                                P["r2as"][i],
                                P["r2bs"][i],
                                es.offset[k : k + CHUNK_SIZE],
                                es.T,
                                es.v1,
                                es.v1err,
                                es.field,
                            )
                        )
            chunk_results = list(self.executor.map(matrix_calc_worker_chunk, tasks))
            results_iter = iter([x for chunk in chunk_results for x in chunk])

        residuals_list = []
        for i, res in active:
            for es in res.estSpecs:
                offsets = np.array(es.offset)
                if use_parallel:
                    est_calc = np.array([next(results_iter) for _ in offsets])
                else:
                    est_calc = self.calc(P, i, offsets, es)
                est_calc = np.nan_to_num(est_calc, nan=1e6, posinf=1e6, neginf=-1e6)
                std_devs = np.array(es.intstd)
                std_devs[std_devs == 0] = 1.0
                residuals_list.extend((np.array(es.int) - est_calc) / std_devs)

        residuals = np.array(residuals_list, dtype=float)
        if len(residuals) == 0:
            if self.verbose:
                print(
                    "Warning (errFunc): No residuals generated. Check active residues and data."
                )
            return (
                np.array([1e6], dtype=float)
                if len(p_flat) > 0
                else np.array([], dtype=float)
            )

        self.chi2 = np.sum(residuals**2)
        self.npar = len(p_flat)
        self.nvar = len(residuals)
        self.dof = max(1, self.nvar - self.npar)
        return residuals

    def fit(self, p0=None, fitting_config=None):
        """
        Performs model fitting. If p0 is None, initial parameters are generated
        from fitting_config via a grid search. Delegates to fit module.
        """
        if fitting_config is not None:
            self.selMethod(fitting_config)

        executor_context = None
        if self.method == "Matrix":
            if self.verbose:
                print("Initializing ProcessPoolExecutor for Matrix method...")
            executor_context = concurrent.futures.ProcessPoolExecutor()
            self.executor = executor_context

        try:
            if p0 is None:
                if fitting_config is None:
                    raise ValueError(
                        "fitting_config must be provided to est_model.fit if p0 is not specified."
                    )
                p0 = fit_module.generate_initial_parameters(self, fitting_config)
            return fit_module.perform_least_squares_fit(self, p0)
        finally:
            if executor_context:
                executor_context.shutdown()
                self.executor = None

    # --- Output: PDFs ---

    def _plot_residues(self, pdfFileName, P=None):
        """One page per active residue: experimental points, plus calc curves if P given."""
        pdf = PdfPages(pdfFileName)
        colorsSet = ["b", "g", "r", "c", "m", "y", "k"]
        for i, res in enumerate(self.dataset.res):
            if not res.active:
                continue
            if not res.estSpecs:
                if self.verbose:
                    print(
                        f"No spectra for active residue {res.label}, skipping PDF page."
                    )
                continue

            fig = figure(figsize=(8.5, 6))
            ax = fig.add_subplot(1, 1, 1)
            all_offsets, all_ints, all_stds = [], [], []

            for j, ep in enumerate(res.estSpecs):
                if not ep.offset:
                    continue
                all_offsets.extend(ep.offset)
                all_ints.extend(ep.int)
                all_stds.extend(ep.intstd)
                c = colorsSet[j % len(colorsSet)]
                ax.errorbar(
                    ep.offset,
                    ep.int,
                    yerr=ep.intstd,
                    fmt=f"{c}o",
                    markersize=3,
                    label=f"exp {ep.v1:.1f}Hz {ep.T * 1000.0:.1f}ms",
                )
                if P is not None:
                    tx = np.array(sorted(ep.offset))
                    ty = np.nan_to_num(np.atleast_1d(self.calc(P, i, tx, ep)))
                    ax.plot(
                        tx,
                        ty,
                        f"{c}-",
                        label=f"calc {ep.v1:.1f}Hz {ep.T * 1000.0:.1f}ms",
                    )

            ax.set_title(res.label)
            ax.set_xlabel("Chemical Shift Offset (ppm)")
            ax.set_ylabel("Intensity")
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width * 0.75, box.height])
            ax.legend(loc="center left", bbox_to_anchor=(1, 0.5), fontsize="small")
            ax.grid(True)

            if all_offsets:
                ax.set_xlim(np.min(all_offsets), np.max(all_offsets))
            if all_ints and all_stds:
                min_y = np.min(np.array(all_ints) - np.array(all_stds))
                max_y = np.max(np.array(all_ints) + np.array(all_stds))
                plot_min_y = min(0.0, min_y * 0.95 if min_y >= 0 else min_y * 1.05)
                plot_max_y = max_y * 1.05 if max_y >= 0 else max_y * 0.95
                if plot_max_y <= plot_min_y:
                    plot_max_y = plot_min_y + 0.1
                ax.set_ylim(plot_min_y, plot_max_y)

            fig.tight_layout(rect=[0, 0, 0.85, 1])
            pdf.savefig(fig)
            close(fig)
        pdf.close()

    def pdf(self, p1_flat, pdfFileName):
        """PDF report of fit results: data + calculated curves per residue."""
        try:
            P = self.seParam(np.asarray(p1_flat, dtype=float))
        except Exception as e:
            print(
                f"Error in seParam during pdf generation: {e}. Skipping PDF generation."
            )
            return
        self._plot_residues(pdfFileName, P)

    def datapdf(self, pdfFileName):
        """PDF of the experimental data only."""
        self._plot_residues(pdfFileName)

    # --- Output: text logs ---

    def _log_header(self, title, extra=()):
        L = ["*" * 51, title, "*" * 51 + "\n"]
        try:
            L.append(f"User: {getlogin()}@{uname()[1]}")
        except Exception:
            L.append("User: Unknown")
        L.append(f"{ctime()}")
        L.extend(extra)
        L.append("*" * 51)
        return L

    def _log_data_tables(self, P, col_hdr="CalcI", width=12):
        """Per-residue tables of experimental vs calculated intensities."""
        L = []
        for ir, ro in enumerate(self.dataset.res):
            if not ro.active:
                continue
            L.append(f"\n# {ro.label}")
            for es in ro.estSpecs:
                L.extend(
                    [
                        f"B0 [MHz]: {es.field:8.3f}",
                        f"T [ms]: {es.T * 1e3:8.3f}",
                        f"v1 [Hz]: {es.v1:8.3f}",
                    ]
                )
                L.append(
                    f"{'Offset':>8} {'ExpI':>12} {'ExpStd':>12} {col_hdr:>{width}}"
                )
                for ko, ovl in enumerate(es.offset):
                    ec = np.nan_to_num(self.calc(P, ir, ovl, es))
                    L.append(
                        f"{ovl:8.3f} {es.int[ko]:12.3f} {es.intstd[ko]:12.3f} {ec:{width}.3f}"
                    )
            L.append("*" * 51)
        return L

    def _log_stats(self, p_ref, mc=False):
        """Chi2 statistics block; recomputes errFunc at p_ref for consistency."""
        _ = self.errFunc(p_ref)
        red = self.chi2 / self.dof if self.dof > 0 else float("inf")
        if mc:
            return [
                "\n\n\n" + "*" * 41,
                f"Results using {self.method} (MC Stats)\n",
                f"Chi2 (init fit): {self.chi2:8.6f}",
                f"red_chi2 (init fit): {red:8.6f}",
                f"dof (init fit):      {self.dof:d}",
                f"nv (init fit):       {self.nvar:d}",
                f"np (init fit):       {self.npar:d}",
                "*" * 43,
                "Fitted Params (Mean from MC +/- StdDev from MC):",
            ]
        return [
            "\n\n\n" + "*" * 41,
            f"Results using {self.method} method\n",
            f"Chi2:     {self.chi2:8.6f}",
            f"red_chi2: {red:8.6f}",
            f"dof:      {self.dof:d}",
            f"nv:       {self.nvar:d}",
            f"np:       {self.npar:d}",
            "*" * 43,
            "Fitted Parameters (+/- StdDev if available):",
        ]

    def _log_params(self, values, stds):
        """Fitted parameter list: global rate constants, then per-residue blocks."""
        gkeys, _, _ = PARAM_LAYOUT[self.method]
        labels = RES_LOG_LABELS[self.method]
        L = []
        pi = 0
        for name in gkeys:
            L.append(f"{name} [s-1]: {values[pi]:8.3f} +/- {stds[pi]:8.3f}")
            pi += 1
        for ro in self.dataset.res:
            if not ro.active:
                continue
            L.extend(["*" * 43, f"residue: {ro.label}"])
            if pi + len(labels) > len(values):
                L.append(f"Error: Parameter index out of bounds for {self.method}.")
                continue
            for lab in labels:
                L.append(f"{lab}{values[pi]:8.4f} +/- {stds[pi]:8.4f}")
                pi += 1
        L.append("*" * 43)
        return L

    def getLogBuffer(self, fit_output_tuple):
        """Text report for a single fit: data tables, chi2 stats, fitted parameters."""
        p_opt, covar = fit_output_tuple[0], fit_output_tuple[1]
        parstd = np.full_like(p_opt, np.nan)
        if (
            covar is not None
            and covar.shape == (len(p_opt), len(p_opt))
            and not np.all(np.isnan(covar))
        ):
            d = np.diag(covar)
            ok = ~np.isnan(d) & ~np.isinf(d) & (d >= 0)
            parstd[ok] = np.sqrt(d[ok])
        elif self.verbose:
            print(
                "Warning: Covar matrix issue or unavailable; std devs are NaN for getLogBuffer."
            )

        P = self.seParam(p_opt)
        L = self._log_header(self.programName)
        L += self._log_data_tables(P)
        L += self._log_stats(p_opt)
        L += self._log_params(p_opt, parstd)
        return "\n".join(L)

    def getLogBufferMC(self, p_center, all_mc_ps_list):
        """Text report for a Monte Carlo run. Returns (log string, mean MC parameters)."""
        if not all_mc_ps_list:
            return "\n".join(
                [
                    "*" * 51,
                    f"{self.programName} - Monte Carlo Analysis - NO SUCCESSFUL RUNS",
                    "*" * 51 + "\n",
                ]
            ), np.array([])

        arr = np.array(all_mc_ps_list)
        mean_p = np.mean(arr, axis=0)
        std_p = np.std(arr, axis=0)

        P = self.seParam(mean_p)
        L = self._log_header(
            f"{self.programName} - Monte Carlo Analysis",
            (f"Num MC runs: {len(all_mc_ps_list)}",),
        )
        L += self._log_data_tables(P, col_hdr="CalcI_MC_Mean", width=18)
        L += self._log_stats(p_center, mc=True)

        gkeys, rkeys, _ = PARAM_LAYOUT[self.method]
        expected = len(gkeys) + sum(r.active for r in self.dataset.res) * len(rkeys)
        if len(mean_p) != expected:
            L.append(
                f"Error: MC parameter array length mismatch. Mean: {len(mean_p)}, StdDev: {len(std_p)}, Expected: {expected}"
            )
            L.append("*" * 43)
        else:
            L += self._log_params(mean_p, std_p)
        return "\n".join(L), mean_p
