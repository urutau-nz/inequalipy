# Inequipy
[Github Link](https://github.com/urutau-nz/inequipy)

Inequipy contains functions for the calculating the inequality of a distribution:
* Kolm-Pollak Equally-Distributed Equivalent (EDE) and Index
* Atkinson EDE and Index
* Gini Index
The Atkinson and Gini approaches are traditionally used for evaluating the inequality of income distribution. The Kolm-Pollak has recently been modified so that it is suitable for quantities that are undesirable, e.g., exposure to health risks or environmental burdens. This means it is suitable for use in urban planning contexts.

When using this code, please cite the following two papers:
* Sheriff, G., & Maguire, K. B. (2020). Health Risk, Inequality Indexes, and Environmental Justice. _Risk Analysis: An Official Publication of the Society for Risk Analysis._ https://doi.org/10.1111/risa.13562
* Logan, T. M., Anderson, M. J., Williams, T., & Conrow, L. (Under review). Measuring inequalities in urban systems: An approach for evaluating the distribution of amenities and burdens. _Computers, Environmental, and Urban Systems_.

### This library contains the following functions:
* `kolm_pollak_ede(a, beta, kappa, weights)` for calculating the Kolm Pollak equally-distributed equivalent
* `kolm_pollak_index(a, beta, kappa, weights)` for calculating the Kolm Pollak inequality index
* `atkinson_ede(a, epsilon, weights)` for calculating the Atkinson equally-distributed equivalent
* `atkinson_index(a, epsilon, weights)` for calculating the Atkinson inequality index
* `gini(a, weights)` for calculating the gini index

### Usage
##### Installation
`pip install inequipy`
##### Usage
Import the package  
`from inequipy import inequipy as ineq`  
Call the required function  
`output = ineq.kolm_pollak_ede(a, beta, kappa, weights)`

### Examples
Check out example.ipynb for examples or https://github.com/MitchellAnderson112/access_inequality_index for the function applied in a real analysis.
