import numpy as np

def ede(a, epsilon = 0.5, weights = None):
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
    sum_atk = 0
    if not weights:
        N = len(a)
    else:
        N = sum(weights)
    x_mean = np.average(a, weights = weights)
    count = 0
    for i in a:
        if not weights:
            sum_atk += i**(1 - epsilon)
        else:
            sum_atk += (i**(1 - epsilon)) * weights[count]
        count += 1
    ede = (sum_atk / N)**(1 / (1 - epsilon))
    return(ede)


def index(a, epsilon = 0.5, weights = None):
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
    ede_atk = ede(a, epsilon, weights)
    x_mean = np.average(a, weights = weights)

    return(1 - (ede_atk / x_mean))
