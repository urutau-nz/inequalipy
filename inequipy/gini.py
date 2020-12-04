import numpy as np

def index(a, weights = None):
    """
    Compute the Gini Coefficient.
    Thanks Gaëtan de Menten: https://stackoverflow.com/a/49571213/5890574

    Parameters
    ----------
    a : array_like
        1-D array containing the values of the distribution.
    weights : array_like, optional
        1-D array of integer weights associated with the values in `a`. Each value in
        `a` contributes to the average according to its associated weight.
        If `weights=None`, then all data in `a` are assumed to have a
        weight equal to one.

    Returns
    -------
    gini : float
        Returns the Gini Coefficient of the distribution provided.
    """
    a = np.asarray(a)
    if weights is not None:
        weights = np.asarray(weights)
        sorted_indices = np.argsort(a)
        sorted_x = a[sorted_indices]
        sorted_w = weights[sorted_indices]
        # Force float dtype to avoid overflows
        cumw = np.cumsum(sorted_w, dtype=float)
        cumxw = np.cumsum(sorted_x * sorted_w, dtype=float)
        return (np.sum(cumxw[1:] * cumw[:-1] - cumxw[:-1] * cumw[1:]) /
                (cumxw[-1] * cumw[-1]))
    else:
        sorted_x = np.sort(a)
        n = len(a)
        cumx = np.cumsum(sorted_x, dtype=float)
        # The above formula, with all weights equal to 1 simplifies to:
        return (n + 1 - 2 * np.sum(cumx) / cumx[-1]) / n
