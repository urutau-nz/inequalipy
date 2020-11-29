import numpy as np
from scipy.integrate import simps

def kolm_pollak_ede(a, beta = None, kappa = None, weights = None):
    """
    Params:
    a: Distribution of data; List
    beta: Set inequality aversion parameter; Integer
    kappa: Inequality aversion parameter; Integer > 0
    weights: For weighted average; List (Length of a)
    
    Returns: Kolm Pollak EDE
    """
    a = np.asanyarray(a)
    if kappa is None:
        if beta is None:
            raise TypeError("you must provide either a beta or kappa aversion parameter")
        kappa = calc_kappa(a, beta, weights)
    if weights is None:
        ede_sum = np.exp(a*-kappa).sum()
        N = len(a)
    else:
        ede_sum = np.multiply(np.exp(a*-kappa), weights).sum()
        N = sum(weights)
    ede = (-1 / kappa) * np.log(ede_sum / N)
    return(ede)


def kolm_pollak_index(a, beta = None, kappa = None, weights = None):
    """
    Params:
    a: Distribution of data; List
    beta: Set inequality aversion parameter; Integer
    kappa: Inequality aversion parameter; Integer > 0
    weights: For weighted average; List (Length of a)
    
    Returns: Kolm Pollak Inequality Index
    """
    if weights is None:
        x_mean = np.mean(a)
    else:
        x_mean = np.average(a, weights = weights)
    a = a - x_mean
    inequality_index = kolm_pollak_ede(a, beta = beta, kappa = kappa, weights = weights)
    return(inequality_index)


def calc_kappa(a, beta, weights = None):
    """
    Params:
        a: Distribution of data; List
        beta: Set inequality aversion parameter; Integer
        weights: For weighted average; List (Length of a)

    Returns: Kappa value calculated by minimising the sum of squares
    """
    if weights is None:
        x_sum = sum(a)
        x_sq_sum = (np.array(a)**2).sum()
    else:
        x_sum = np.multiply(a, weights).sum()
        x_sq_sum = np.multiply(a**2, weights).sum()
    kappa = beta * (x_sum / x_sq_sum)
    return(kappa)


####
# Other inequality functions
####

def atkinson_ede(a, epsilon = -0.5, weights = None):
    """
    Params:
        a: Distribution of data; List
        epsilon: Set inequality aversion parameter; Integer
        weights: For weighted average; List (Length of a)

    Returns: Normal Atkinson EDE

    NOTE: this is meant for distributions where high values are better
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

def atkinson_index(a, epsilon = -0.5, weights = None):
    """
    Params:
        a: Distribution of data; List
        epsilon: Set inequality aversion parameter; Integer
        weights: For weighted average; List (Length of a)

    Returns: Atkinson index
    """
    ede_sum = 0 # init sum
    if not weights:
        N = len(a)
    else:
        N = sum(weights) #sum of data
    x_mean = np.average(a, weights = weights)
    count = 0

    for i in a:
        if not weights:
            ede_sum += (i**(1-epsilon))
        else:
            ede_sum += (i**(1-epsilon)) * weights[count]
        count += 1

    ede = (ede_sum / N)**(1 / (1 - epsilon))

    index = 1 - (ede / x_mean)
    return(index)

def gini(a, beta = -0.5, weights = None):
    """
    Params:
        a: Distribution of data; List
        beta: Set inequality aversion parameter; Integer
        weights: For weighted average; List (Length of a)
        
    Returns: Gini index
    """
    # based on bottom eq:
    # http://www.statsdirect.com/help/generatedimages/equations/equation154.svg
    # from:
    # http://www.statsdirect.com/help/default.htm#nonparametric_methods/gini.htm
    # All values are treated equally, arrays must be 1d:
    if weights:
        array = np.repeat(array,weights)
    array = array.flatten()
    if np.amin(array) < 0:
        # Values cannot be negative:
        array -= np.amin(array)
    # Values cannot be 0:
    array += 0.0000001
    # Values must be sorted:
    array = np.sort(array)
    # Index per array element:
    index = np.arange(1,array.shape[0]+1)
    # Number of array elements:
    n = array.shape[0]
    # Gini coefficient:
    return ((np.sum((2 * index - n  - 1) * array)) / (n * np.sum(array)))

