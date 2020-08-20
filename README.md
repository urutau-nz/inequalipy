# Inequality Metrics
[Github Link](https://github.com/michael-j-freeman/urutau-inequality-metrics)

Inequality Metrics contains functions for the Kolm-Pollak, Atkinson EDE and Gini Index aproaches to calculating inequality distributions. These methods allow for better statistical modeling of inequality than other, more standard methods as they include a defined "Inequality Aversion Parameter" which allows for both positive and negative distributions to be easily calculated.
<br/>

### This library contains the following functions:<br/>
* `kolm_pollak_ede()` for calculating the Kolm Pollak equally distributed equivilant
* `kolm_pollak_index()` for calculating the Kolm Pollak inequality index
* `atkinson_ede()` for calcuating the Atkinson equally distributed equivilant
* `atkinson_index()` for calcuating the Atkinson inequality index
* `gini_index()` for calculating the gini index

### Usage

##### Installation
`pip install InequalityMetrics`
##### Usage
Import the package
`import InequalityMetrics`
Call the required function
`output = InequalityMetrics.kolm_pollak_ede()`