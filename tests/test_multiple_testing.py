import numpy as np

from pangwasflow.multiple_testing import benjamini_hochberg


def test_bh_bounds_and_monotonicity():
    pvalues = np.array([0.001, 0.02, 0.03, 0.5, 0.9])
    adjusted = benjamini_hochberg(pvalues)
    assert np.all((adjusted >= 0) & (adjusted <= 1))
    assert np.all(np.diff(adjusted) >= 0)
