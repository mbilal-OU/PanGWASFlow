from __future__ import annotations

import numpy as np


def benjamini_hochberg(pvalues: np.ndarray) -> np.ndarray:
    """Return Benjamini-Hochberg adjusted P values."""
    pvalues = np.asarray(pvalues, dtype=float)
    output = np.full_like(pvalues, np.nan, dtype=float)
    finite = np.isfinite(pvalues)
    observed = pvalues[finite]
    if observed.size == 0:
        return output

    order = np.argsort(observed)
    ranked = observed[order]
    adjusted = ranked * observed.size / np.arange(1, observed.size + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0.0, 1.0)

    restored = np.empty_like(adjusted)
    restored[order] = adjusted
    output[finite] = restored
    return output
