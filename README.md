# Inequalipy
[Github Link](https://github.com/urutau-nz/inequalipy)  
[Python Pip Link](https://pypi.org/project/inequalipy)  
[GitHub Link for R version](https://github.com/urutau-nz/kolmpollak-R)

Inequalipy contains functions for the calculating the inequality of a distribution:
* Kolm-Pollak Equally-Distributed Equivalent (EDE) and Index
* Atkinson EDE and Index
* Gini Index  

The Atkinson and Gini approaches are traditionally used for evaluating the inequality of income distribution. The Kolm-Pollak has recently been modified so that it is suitable for quantities that are undesirable, e.g., exposure to health risks or environmental burdens. This means it is suitable for use in urban planning contexts.  
The Kolm-Pollak functions are also separable, meaning they can be used to evaluate
subgroups and thus evaluate the inequity of the distribution.

When using this code, please cite the following two papers:
* Sheriff, G., & Maguire, K. B. (2020). Health Risk, Inequality Indexes, and Environmental Justice. _Risk Analysis: An Official Publication of the Society for Risk Analysis._ https://doi.org/10.1111/risa.13562
* Logan, T. M., Anderson, M. J., Williams, T. G., & Conrow, L. (2021). Measuring inequalities in urban systems: An approach for evaluating the distribution of amenities and burdens. Computers, Environment and Urban Systems, 86, 101590. https://doi.org/10.1016/j.compenvurbsys.2020.101590

### This library contains the following functions:
* `kolmpollak.ede(a, epsilon, kappa, weights)` for calculating the Kolm-Pollak equally-distributed equivalent (EDE)
* `kolmpollak.index(a, epsilon, kappa, weights)` for calculating the Kolm-Pollak inequality index
* `atkinson.ede(a, epsilon, weights)` for calculating the Atkinson equally-distributed equivalent
* `atkinson.index(a, epsilon, weights)` for calculating the Atkinson inequality index
* `gini(a, weights)` for calculating the Gini index

### Usage
##### Installation
`pip install inequalipy`
##### Usage
Import the package and call the required function:
```
import inequalipy as ineq  
ineq.kolmpollak.ede(a, epsilon, weights)
ineq.gini(a)
```  
Alternatively, you may import all of the packages:
```
from inequalipy import *
gini(a)
atkinson.index(a)
```

### Examples
Check out example.ipynb for examples or https://github.com/MitchellAnderson112/access_inequality_index for the function applied in a non-trivial context.
