import numpy as np

def ede(a, epsilon = None, kappa = None, weights = None):
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
    return(-1 / kappa) * np.log(ede_sum / N)


def index(a, epsilon = None, kappa = None, weights = None):
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

    return ede(a - x_mean, epsilon = epsilon, kappa = kappa, weights = weights)


def calc_kappa(a, epsilon, weights = None):
    """
    Converts the inequality aversion parameter used in Atkinson's formulae (epsilon)
    into the form for the Kolm-Pollak formulae (kappa).
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
