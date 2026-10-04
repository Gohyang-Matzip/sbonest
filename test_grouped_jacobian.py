"""Regression checks for grouped finite differences with fixed parameters."""

import numpy as np
from scipy.optimize._numdiff import approx_derivative

from fit import _block_jacobian


SIZES = [8, 11, 9]
N_GLOBAL, N_LOCAL = 3, 7
PARAMETERS = np.linspace(-1.2, 2.4, N_GLOBAL + len(SIZES) * N_LOCAL)


def residual(full):
    """Nonlinear blocks share globals but have independent residue locals."""
    blocks = []
    globals_ = full[:N_GLOBAL]
    for residue, size in enumerate(SIZES):
        local = full[N_GLOBAL + residue * N_LOCAL:N_GLOBAL + (residue + 1) * N_LOCAL]
        grid = np.linspace(0.15, 0.95, size)
        block = sum(
            (key + 1) * np.sin(grid * value) + 0.1 * value**2
            for key, value in enumerate(globals_)
        )
        block += sum(
            np.cos((key + 1) * grid) * (value + 0.1 * value**2)
            for key, value in enumerate(local)
        )
        block += globals_[0] * local[0] * grid
        blocks.append(block)
    return np.concatenate(blocks)


def reduced_residual(free, initial=PARAMETERS):
    def evaluate(values):
        full = initial.copy()
        full[free] = values
        return residual(full)

    return evaluate


def check_legacy_derivative():
    """No new arguments retains the original one-sided result and call count."""
    calls = []

    def counted(values):
        calls.append(values.copy())
        return residual(values)

    actual = _block_jacobian(counted, SIZES, N_GLOBAL, N_LOCAL)(PARAMETERS)
    assert len(calls) == 11, len(calls)
    base = residual(PARAMETERS)
    reference = np.empty_like(actual)
    for column, value in enumerate(PARAMETERS):
        shifted = PARAMETERS.copy()
        shifted[column] += np.sqrt(np.finfo(float).eps) * np.sign(value) * max(1, abs(value))
        reference[:, column] = (residual(shifted) - base) / (shifted[column] - value)
    np.testing.assert_array_equal(actual, reference)


def check_free_order_and_inactive_groups():
    subsets = [
        np.arange(len(PARAMETERS)),
        np.delete(np.arange(len(PARAMETERS)), N_GLOBAL + 1),
        np.array([23, 0, 5, 16, 2, 11, 3]),
        np.array([23, 5, 3]),  # No globals and no free locals for residue 1.
        np.array([2, 0]),  # Every local group is inactive.
    ]
    for free in subsets:
        evaluate = reduced_residual(free)
        calls = []

        def counted(values):
            calls.append(values.copy())
            return evaluate(values)

        values = PARAMETERS[free]
        actual = _block_jacobian(
            counted, SIZES, N_GLOBAL, N_LOCAL, relative_step=1e-6, free=free
        )(values)
        expected_calls = 1 + np.count_nonzero(free < N_GLOBAL)
        expected_calls += len(set((free[free >= N_GLOBAL] - N_GLOBAL) % N_LOCAL))
        assert len(calls) == expected_calls <= 11, (free, len(calls), expected_calls)
        assert actual.shape == (sum(SIZES), len(free)), actual.shape
        reference = approx_derivative(evaluate, values, method="3-point")
        np.testing.assert_allclose(actual, reference, rtol=2e-5, atol=2e-6)


def check_bounded_derivative():
    # Reordered globals/locals include negative lower bounds, positive upper
    # bounds, near bounds, and a narrow interval requiring a shorter step.
    free = np.array([23, 0, 5, 16, 2, 11, 3])
    values = np.array([1.0, -0.5, 0.6, -0.4, 0.25, 0.2, -0.3])
    lower = np.array([0.0, -0.5, 0.0, -0.4 - 1e-9, 0.25 - 2e-7, -1.0, -1.0])
    upper = np.array([1.0, 1.0, 0.6 + 1e-9, 1.0, 0.25 + 1e-7, 0.2, 1.0])
    initial = PARAMETERS.copy()
    initial[free] = values
    evaluate = reduced_residual(free, initial)
    calls = []

    def bounded(values):
        assert np.all(values >= lower) and np.all(values <= upper), values
        calls.append(values.copy())
        return evaluate(values)

    actual = _block_jacobian(
        bounded, SIZES, N_GLOBAL, N_LOCAL,
        relative_step=1e-6, free=free, bounds=(lower, upper),
    )(values)
    reference = approx_derivative(evaluate, values, method="3-point", bounds=(lower, upper))
    np.testing.assert_allclose(actual, reference, rtol=2e-5, atol=2e-6)
    assert len(calls) <= 11, len(calls)
    displacements = np.asarray(calls[1:]) - values
    assert np.any(displacements[:, 0] < 0)  # Exact positive upper bound.
    assert np.any(displacements[:, 1] > 0)  # Exact negative lower bound.
    assert np.any(displacements[:, 2] < 0)  # Near positive upper bound.
    assert np.any(displacements[:, 3] > 0)  # Near negative lower bound.
    assert np.any(displacements[:, 4] < 0)  # Shorten toward farther bound.


def check_scalar_bounds_without_free():
    values = np.linspace(0.0, 1.0, len(PARAMETERS))

    def bounded(values):
        assert np.all(values >= 0.0) and np.all(values <= 1.0), values
        return residual(values)

    actual = _block_jacobian(
        bounded, SIZES, N_GLOBAL, N_LOCAL, relative_step=1e-6, bounds=(0.0, 1.0)
    )(values)
    reference = approx_derivative(bounded, values, method="3-point", bounds=(0.0, 1.0))
    np.testing.assert_allclose(actual, reference, rtol=2e-5, atol=2e-6)


def check_three_point_subsets():
    for index, free in enumerate([
        np.delete(np.arange(len(PARAMETERS)), N_GLOBAL + 1),
        np.array([23, 0, 5, 16, 2, 11, 3]),
        np.array([23, 5, 3]),
        np.array([2, 0]),
    ]):
        values = PARAMETERS[free].copy()
        values[-1] = 0.0  # Uses the automatic step even with rel_step specified.
        # Cover both representable tiny steps and relative-step underflow.
        values[0] = 1e-310 if index % 2 == 0 else 1e-320
        initial = PARAMETERS.copy()
        initial[free] = values
        evaluate = reduced_residual(free, initial)
        calls = []

        def counted(values):
            calls.append(values.copy())
            return evaluate(values)

        actual = _block_jacobian(
            counted, SIZES, N_GLOBAL, N_LOCAL, relative_step=1e-5,
            free=free, method="3-point",
        )(values)
        reference = approx_derivative(evaluate, values, method="3-point", rel_step=1e-5)
        np.testing.assert_array_equal(actual, reference)
        group_count = np.count_nonzero(free < N_GLOBAL)
        group_count += len(set((free[free >= N_GLOBAL] - N_GLOBAL) % N_LOCAL))
        assert len(calls) == 1 + 2 * group_count <= 21, len(calls)


def check_three_point_bounds_and_steps():
    # Columns 23, 16 and 9 share one local group but need backward, forward
    # and central stencils, respectively. Other columns exercise tightened
    # central stencils, a zero fallback, and a relative step below unit scale.
    free = np.array([23, 0, 9, 16, 2, 11, 3])
    values = np.array([1.0, 0.25, 0.6, -0.4, 0.0, 0.2, -0.3])
    lower = np.array([0.0, 0.25 - 2e-7, 0.0, -0.4, 0.0, -1.0, -1.0])
    upper = np.array([1.0, 0.25 + 1e-7, 1.0, 1.0, 1.0, 0.2, 1.0])
    initial = PARAMETERS.copy()
    initial[free] = values
    evaluate = reduced_residual(free, initial)
    for relative_step in [1e-5, None, 0.0, 1e-20]:
        calls = []

        def bounded(values):
            assert np.all(values >= lower) and np.all(values <= upper), values
            calls.append(values.copy())
            return evaluate(values)

        actual = _block_jacobian(
            bounded, SIZES, N_GLOBAL, N_LOCAL, relative_step=relative_step,
            free=free, bounds=(lower, upper), method="3-point",
        )(values)
        reference_calls = []

        def reference_residual(values):
            reference_calls.append(values.copy())
            return evaluate(values)

        reference = approx_derivative(
            reference_residual, values, method="3-point", rel_step=relative_step,
            bounds=(lower, upper),
        )
        # The independent dense stencil can produce roundoff on exactly
        # constant residual blocks; the grouped result retains structural zero.
        for column, original in enumerate(free):
            if original < N_GLOBAL:
                affected = slice(None)
            else:
                residue = (original - N_GLOBAL) // N_LOCAL
                start = sum(SIZES[:residue])
                affected = slice(start, start + SIZES[residue])
                unaffected = np.ones(sum(SIZES), dtype=bool)
                unaffected[affected] = False
                np.testing.assert_array_equal(actual[unaffected, column], 0.0)
            np.testing.assert_array_equal(actual[affected, column], reference[affected, column])
            actual_positions = sorted(v[column] for v in calls[1:] if v[column] != values[column])
            expected_positions = sorted(v[column] for v in reference_calls[1:] if v[column] != values[column])
            np.testing.assert_array_equal(actual_positions, expected_positions)


def check_mixed_local_stencils():
    sizes = [2, 1, 3, 1, 2, 1, 2, 1]
    values = np.array([0.3, -0.4, -0.5, -0.6, 0.25, 0.25, -0.25, 0.0])
    limits = np.array([
        [-1.0, 1.0], [-1.0, 1.0], [-0.5, 1.0], [-1.0, -0.6],
        [0.25 - 2**-24, 0.25 + 2**-25],
        [0.25 - 2**-24, 0.25 + 2**-28],
        [-0.25 - 2**-28, -0.25 + 2**-24], [-1.0, 1.0],
    ])
    lower, upper = limits.T
    calls = []

    def evaluate(values):
        assert np.all(values >= lower) and np.all(values <= upper), values
        calls.append(values.copy())
        return np.concatenate([
            np.arange(1, size + 1) * value**3 + 0.5 * value**2 + value
            for value, size in zip(values, sizes)
        ])

    actual = _block_jacobian(
        evaluate, sizes, 0, 1, relative_step=1e-5,
        bounds=(lower, upper), method="3-point",
    )(values)
    assert len(calls) == 3, len(calls)
    grouped_calls = calls.copy()
    reference = approx_derivative(
        evaluate, values, method="3-point", rel_step=1e-5, bounds=(lower, upper)
    )
    for column, size in enumerate(sizes):
        start = sum(sizes[:column])
        affected = slice(start, start + size)
        np.testing.assert_array_equal(actual[affected, column], reference[affected, column])
        # Dense reference calls are base followed by two evaluations per column.
        reference_positions = [calls[4 + 2 * column][column], calls[5 + 2 * column][column]]
        np.testing.assert_array_equal(
            [grouped_calls[1][column], grouped_calls[2][column]], reference_positions
        )


if __name__ == "__main__":
    for check in [
        check_legacy_derivative,
        check_free_order_and_inactive_groups,
        check_bounded_derivative,
        check_scalar_bounds_without_free,
        check_three_point_subsets,
        check_three_point_bounds_and_steps,
        check_mixed_local_stencils,
    ]:
        check()
        print("PASS", check.__name__, flush=True)
