# Inequality Metrics
[Github Link](https://github.com/urutau-nz/inequality-metrics)

Inequality Metrics contains functions for the Kolm-Pollak, Atkinson EDE and Gini Index aproaches to calculating inequality distributions. These methods allow for better statistical modeling of inequality than other, more standard methods as they include a defined "Inequality Aversion Parameter" which allows for both positive and negative distributions to be easily calculated.
<br/>

### This library contains the following functions:<br/>
* `kolm_pollak_ede(a, beta, kappa, weights)` for calculating the Kolm Pollak equally distributed equivilant
    *   Parameters:
        * a: Distribution of data; List
        * beta: Set inequality aversion parameter; Integer|Optional 
        * kappa: Inequality aversion parameter; Integer > 0|Optional 
        * weights: For weighted average; List (Length of a)|Optional 
* `kolm_pollak_index(a, beta, kappa, weights)` for calculating the Kolm Pollak inequality index
    *   Parameters:
        * a: Distribution of data; List
        * beta: Set inequality aversion parameter; Integer|Optional 
        * kappa: Inequality aversion parameter; Integer > 0|Optional 
        * weights: For weighted average; List (Length of a)|Optional 
* `atkinson_ede(a, epsilon, weights)` for calcuating the Atkinson equally distributed equivilant
    *   Parameters:
        *    a: Distribution of data; List
        *   epsilon: Set inequality aversion parameter; Integer|Optional 
        *   weights: For weighted average; List (Length of a)|Optional 
* `atkinson_index(a, epsilon, weights)` for calcuating the Atkinson inequality index
     *   Parameters:
         *    a: Distribution of data; List
         *   epsilon: Set inequality aversion parameter; Integer|Optional 
         *   weights: For weighted average; List (Length of a)|Optional 
* `gini(a, beta, weights)` for calculating the gini index
     *   Parameters:
         *    a: Distribution of data; List
         *    beta: Set inequality aversion parameter; Integer|Optional 
         *    weights: For weighted average; List (Length of a)|Optional 

### Usage
##### Installation
`pip install InequalityMetrics`
##### Usage
Import the package </br>
`from InequalityMetrics import InequalityMetrics`<br/>
Call the required function<br/>
`output = InequalityMetrics.kolm_pollak_ede(a, beta, kappa, weights)`<br/>