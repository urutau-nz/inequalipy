import numpy as np
from scipy.integrate import simps

def kolm_pollak_ede(a, epsilon = None, kappa = None, weights = None):
    """
    Compute the Kolm-Pollak Equally-Distributed Equivalent (EDE).
    The Kolm-Pollak EDE and Index are suitable for distributions of desirable and
    undesirable quantities. That is, a desirable quantity (like income) is where
    having more of the quantity is desirable; compared with an undesirable
    quantity like health risk, where less is better.

    Parameters
    ----------
    a : array_like
        1-D array containing the values of the distribution.
    epsilon : float
        The inequality aversion parameter from the Atkinson formulae.
            If epsilon > 0 then the quantity is desirable (more is better).
    kappa : float
        The inequality aversion parameter from the Kolm-Pollak formulae.
            If kappa > 0 then the quantity is desirable (more is better).
    weights : array_like, optional
        1-D array of integer weights associated with the values in `a`. Each value in
        `a` contributes to the average according to its associated weight.
        If `weights=None`, then all data in `a` are assumed to have a
        weight equal to one.

    Returns
    -------
    value : float
        Returns the Kolm-Pollak Equally-Distributed Equivalent of the distribution provided.
    """
    a = np.asanyarray(a)
    if kappa is None:
        if epsilon is None:
            raise TypeError("you must provide either a epsilon or kappa aversion parameter")
        kappa = calc_kappa(a, epsilon, weights)
    if weights is None:
        ede_sum = np.exp(a*-kappa).sum()
        N = len(a)
    else:
        ede_sum = np.multiply(np.exp(a*-kappa), weights).sum()
        N = sum(weights)
    ede = (-1 / kappa) * np.log(ede_sum / N)
    return(ede)


def kolm_pollak_index(a, epsilon = None, kappa = None, weights = None):
    """
    Compute the Kolm-Pollak Index.
    The Kolm-Pollak EDE and Index are suitable for distributions of desirable and
    undesirable quantities. That is, a desirable quantity (like income) is where
    having more of the quantity is desirable; compared with an undesirable
    quantity like health risk, where less is better.

    Parameters
    ----------
    a : array_like
        1-D array containing the values of the distribution.
    epsilon : float
        The inequality aversion parameter from the Atkinson formulae.
            If epsilon > 0 then the quantity is desirable (more is better).
    kappa : float
        The inequality aversion parameter from the Kolm-Pollak formulae.
            If kappa > 0 then the quantity is desirable (more is better).
    weights : array_like, optional
        1-D array of integer weights associated with the values in `a`. Each value in
        `a` contributes to the average according to its associated weight.
        If `weights=None`, then all data in `a` are assumed to have a
        weight equal to one.

    Returns
    -------
    value : float
        Returns the Kolm-Pollak Index of the distribution provided.
    """
    if weights is None:
        x_mean = np.mean(a)
    else:
        x_mean = np.average(a, weights = weights)

    return kolm_pollak_ede(a - x_mean, epsilon = epsilon, kappa = kappa, weights = weights)


def calc_kappa(a, epsilon, weights = None):
    """
    Converts the inequality aversion parameter used in Atkinson's formulae (epsilon)
    into the form for the Kolm Pollak formulae (kappa).
    If epsilon > 0 then the quantity is desirable (more is better).

    Parameters
    ----------
    a : array_like
        1-D array containing the values of the distribution.
    epsilon : float
        The inequality aversion parameter from the Atkinson formulae.
            If epsilon > 0 then the quantity is desirable (more is better).
    weights : array_like, optional
        1-D array of integer weights associated with the values in `a`. Each value in
        `a` contributes to the average according to its associated weight.
        If `weights=None`, then all data in `a` are assumed to have a
        weight equal to one.

    Returns
    -------
    value : float
        Returns the inequality aversion parameter for the Kolm-Pollak formulae
    """
    if weights is None:
        x_sum = sum(a)
        x_sq_sum = (np.array(a)**2).sum()
    else:
        x_sum = np.multiply(a, weights).sum()
        x_sq_sum = np.multiply(a**2, weights).sum()
    return(epsilon * (x_sum / x_sq_sum))


####
# Other inequality functions
####

def atkinson_ede(a, epsilon = 0.5, weights = None):
    """
    Compute the Atkinson Equally-Distributed Equivalent.
    The Atkinson EDE and Index are only suitable for distributions of desirable
    quantities (where having more of the quantity is desirable), e.g., income.

    Parameters
    ----------
    a : array_like
        1-D array containing the values of the distribution.
    epsilon : float
        The inequality aversion parameter. epsilon > 0.
    weights : array_like, optional
        1-D array of integer weights associated with the values in `a`. Each value in
        `a` contributes to the average according to its associated weight.
        If `weights=None`, then all data in `a` are assumed to have a
        weight equal to one.

    Returns
    -------
    value : float
        Returns the Atkinson EDE of the distribution provided.
    """
    at_sum = 0
    if not weights:
        N = len(a)
    else:
        N = sum(weights)
    x_mean = np.average(a, weights = weights)
    count = 0
    for i in a:
        if not weights:
            at_sum += i**(1 - epsilon)
        else:
            at_sum += (i**(1 - epsilon)) * weights[count]
        count += 1
    ede = (at_sum / N)**(1 / (1 - epsilon))
    return(ede)

def atkinson_index(a, epsilon = 0.5, weights = None):
    """
    Compute the Atkinson Index.
    The Atkinson EDE and Index are only suitable for distributions of desirable
    quantities (where having more of the quantity is desirable), e.g., income.

    Parameters
    ----------
    a : array_like
        1-D array containing the values of the distribution.
    epsilon : float
        The inequality aversion parameter. epsilon > 0.
    weights : array_like, optional
        1-D array of integer weights associated with the values in `a`. Each value in
        `a` contributes to the average according to its associated weight.
        If `weights=None`, then all data in `a` are assumed to have a
        weight equal to one.

    Returns
    -------
    value : float
        Returns the Atkinson Index of the distribution provided.
    """
    ede = atkinson_ede(a, epsilon, weights)
    x_mean = np.average(a, weights = weights)

    return(1 - (ede / x_mean))

def gini(a, weights = None):
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
