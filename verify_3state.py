import numpy as np
import sys
import os

# Add current directory to path
sys.path.append(os.getcwd())

from estmodel import est_model
from est_data import EstDataSet, Residue, EstSpec


def create_mock_dataset():
    ds = EstDataSet()
    res = Residue()
    res.label = "TEST1"
    res.active = True

    spec = EstSpec()
    spec.field = 600.0  # MHz
    spec.T = 0.4  # 400 ms
    spec.v1 = 25.0  # Hz
    spec.v1err = 0.0
    spec.offset = np.linspace(-5, 5, 21).tolist()
    spec.int = np.ones(21).tolist()  # Dummy intensity
    spec.intstd = np.ones(21).tolist()
    spec.initdw = 2.0
    spec.initr2a = 10.0
    spec.initr2b = 20.0

    res.estSpecs.append(spec)
    ds.res.append(res)
    return ds


def test_3state_linear():
    print("\n--- Testing 3-State Linear ---")
    model = est_model()
    model.dataset = create_mock_dataset()
    model.method = "Matrix_3st_Linear"

    # kab, kba, kbc, kcb, dG, dw, dwC, R1, R2a, R2b, R2c
    # Let's simulate a dip at 0 (A), 2 (B), and 4 (C)
    # dG=0, dw=2 (B at 2), dwC=4 (C at 4)

    kab = 20.0
    kba = 20.0
    kbc = 20.0
    kcb = 20.0
    dG = 0.0
    dw = 2.0
    dwC = 4.0
    R1 = 1.0
    R2a = 10.0
    R2b = 20.0
    R2c = 20.0

    offsets = np.linspace(-2, 6, 9)  # -2, -1, 0, 1, 2, 3, 4, 5, 6
    print(f"Offsets: {offsets}")

    intensities = []
    for off in offsets:
        val = model.matrix_calc_3st_linear(
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
            off,
            0.4,
            25.0,
            0.0,
            600.0,
            observable="A",
        )
        intensities.append(val)

    print("Calculated Intensities (Linear):")
    print(np.array(intensities))

    # Check if we have dips at 0, 2, 4
    # Intensity should be lower at these offsets
    # 0 is index 2, 2 is index 4, 4 is index 6

    i_arr = np.array(intensities)
    print(f"I at 0ppm: {i_arr[2]:.4f}")
    print(f"I at 2ppm: {i_arr[4]:.4f}")
    print(f"I at 4ppm: {i_arr[6]:.4f}")
    print(f"I at 6ppm (far): {i_arr[8]:.4f}")

    assert np.isfinite(i_arr).all()
    assert i_arr[2] < i_arr[8] and i_arr[4] < i_arr[8] and i_arr[6] < i_arr[8]
    print("SUCCESS: Dips observed at expected positions.")


def test_3state_triangle():
    print("\n--- Testing 3-State Triangle ---")
    model = est_model()
    model.dataset = create_mock_dataset()
    model.method = "Matrix_3st_Triangle"

    # kab, kba, kbc, kcb, kca, kac
    kab = 20.0
    kba = 20.0
    kbc = 20.0
    kcb = 20.0
    kca = 20.0
    kac = 20.0

    dG = 0.0
    dw = 2.0
    dwC = -2.0  # C at -2
    R1 = 1.0
    R2a = 10.0
    R2b = 20.0
    R2c = 20.0

    offsets = np.linspace(-4, 4, 9)  # -4, -3, -2, -1, 0, 1, 2, 3, 4
    print(f"Offsets: {offsets}")

    intensities = []
    for off in offsets:
        val = model.matrix_calc_3st_triangle(
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
            off,
            0.4,
            25.0,
            0.0,
            600.0,
        )
        intensities.append(val)

    print("Calculated Intensities (Triangle):")
    print(np.array(intensities))

    # Dips at 0 (A), 2 (B), -2 (C)
    # 0 is index 4, 2 is index 6, -2 is index 2
    i_arr = np.array(intensities)
    print(f"I at -2ppm: {i_arr[2]:.4f}")
    print(f"I at 0ppm: {i_arr[4]:.4f}")
    print(f"I at 2ppm: {i_arr[6]:.4f}")
    print(f"I at 4ppm (far): {i_arr[8]:.4f}")

    assert np.isfinite(i_arr).all()
    assert i_arr[2] < i_arr[8] and i_arr[4] < i_arr[8] and i_arr[6] < i_arr[8]
    print("SUCCESS: Dips observed at expected positions.")


def test_fit_pipeline():
    print("\n--- Testing Fit Pipeline (Synthetic Profiles) ---")
    for method, rates in [
        ("Matrix_3st_Linear", [20.0, 40.0, 10.0, 30.0]),
        ("Matrix_3st_Triangle", [20.0, 40.0, 10.0, 30.0, 5.0, 15.0]),
    ]:
        model = est_model()
        model.dataset = create_mock_dataset()
        model.method = method
        spec = model.dataset.res[0].estSpecs[0]
        spec.offset = np.linspace(-4, 6, 41).tolist()
        spec.intstd = [0.01] * 41
        truth = np.array(rates + [0.0, 2.0, 4.0, 1.0, 10.0, 20.0, 20.0])
        spec.int = model.calc(
            model.seParam(truth), 0, np.array(spec.offset), spec
        ).tolist()
        optimized, _ = model.fit(fitting_config={"Method": method})
        chi2 = np.sum(model.errFunc(optimized) ** 2)
        assert np.isfinite(optimized).all() and chi2 < 1e-8, (method, chi2)
        print(f"SUCCESS: {method} synthetic fit, chi2={chi2:.3g}")


if __name__ == "__main__":
    test_3state_linear()
    test_3state_triangle()
    test_fit_pipeline()
