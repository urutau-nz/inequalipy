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

def atkinson_adjusted_ede(a, beta = -0.5, weight = None):
    """
    Params:
        a: Distribution of data; List
        beta: Set inequality aversion parameter; Integer
        weights: For weighted average; List (Length of a)

    Returns: Adjusted Atkinson EDE
    """
    ede_sum = 0
    if not weight:
        N = len(a)
    else:
        N = sum(weight)
    x_mean = np.average(a, weights = weight)
    count = 0
    for i in a:
        if not weight:
            ede_sum += (i**(1-beta))
        else:
            ede_sum += (i**(1-beta)) * weight[count]
        count += 1
    ede = (ede_sum / N)**(1 / (1 - beta))
    return(ede)

def atkinson_adjusted_index(a, beta = -0.5, weight = None):
    """
    Params:
        a: Distribution of data; List
        beta: Set inequality aversion parameter; Integer
        weights: For weighted average; List (Length of a)

    Returns: Adjusted Atkinson index
    """
    ede_sum = 0
    if not weight:
        N = len(a)
    else:
        N = sum(weight)
    x_mean = np.average(a, weights = weight)
    count = 0
    for i in a:
        if not weight:
            ede_sum += (i**(1-beta))
        else:
            ede_sum += (i**(1-beta)) * weight[count]
        count += 1
    ede = (ede_sum / N)**(1 / (1 - beta))
    index = (ede / x_mean) - 1
    return(index)

def atkinson_ede(a, epsilon = -0.5, weight = None):
    """
    Params:
        a: Distribution of data; List
        epsilon: Set inequality aversion parameter; Integer
        weights: For weighted average; List (Length of a)

    Returns: Normal Atkinson EDE

    NOTE: this is meant for distributions where high values are better
    """
    at_sum = 0
    if not weight:
        N = len(a)
    else:
        N = sum(weight)
    x_mean = np.average(a, weights = weight)
    count = 0
    for i in a:
        if not weight:
            at_sum += i**(1 - epsilon)
        else:
            at_sum += (i**(1 - epsilon)) * weight[count]
        count += 1
    ede = (at_sum / N)**(1 / (1 - epsilon))
    return(ede)

def atkinson_index(a, epsilon = -0.5, weight = None):
    """
    Params:
        a: Distribution of data; List
        epsilon: Set inequality aversion parameter; Integer
        weights: For weighted average; List (Length of a)

    Returns: Atkinson index
    """
    ede_sum = 0
    if not weight:
        N = len(a)
    else:
        N = sum(weight)
    x_mean = np.average(a, weights = weight)
    count = 0
    for i in a:
        if not weight:
            ede_sum += (i**(1-epsilon))
        else:
            ede_sum += (i**(1-epsilon)) * weight[count]
        count += 1
    ede = (ede_sum / N)**(1 / (1 - epsilon))
    index = 1 - (ede / x_mean)
    return(index)

def gini_index(a, beta = -0.5, weight = None):
    """
    Params:
        a: Distribution of data; List
        beta: Set inequality aversion parameter; Integer
        weights: For weighted average; List (Length of a)
        
    Returns: Gini index
    """
    area_total = simps(np.arange(0,101,1), dx=1)
    if not weight:
        N = len(a)
    else:
        N = sum(weight)
    a = list(np.sort(a))
    a_percent = []
    weight_perc = []
    w_perc_sum = 0
    perc_sum = 0
    if not weight:
        for i in a:
            perc_sum += i/data_tot*100
            a_percent.append(perc_sum)
        area_real = simps(a_percent)
        area_diff = area_total - area_real
        gini = round((area_diff / area_total), 3)
    else:
        weight_tot = sum(weight)
        data_tot = sum(a)
        for i in a:
            perc_sum += i/data_tot*100
            a_percent.append(perc_sum)
        for i in weight:
            w_perc_sum += i/weight_tot*100
            weight_perc.append(w_perc_sum)
        area_real = simps(a_percent, weight_perc)
        area_diff = area_total - area_real
        gini = round((area_diff / area_total), 3)
    return(gini)