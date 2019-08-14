# MaOEA

<!-- badges: start -->
<!-- badges: end -->

The goal of MaOEA is to facilitate easy hybridization of algorithms for many objective optimization. In the package, several algorithms are available: SMS-EMOA, NSGA-III, and MO-CMA-ES. Each of these algorithms can be accessed independently. Using the main function, the algorithms can be called for specific number of iterations. Alternatively, if the hybridization follows a more complex rule, users may prefer to call the algorithm directly in their optimization loop. This will call the algorithm (i.e., the offspring generation and selection scheme) for a single iteration.

## Installation

MaOEA requires the users to have installed Python (see https://www.python.org) and being able to use the PyGMO module. Installation instruction for PyGMO is available in https://esa.github.io/pagmo2/install.html.

After Python installed and PyGMO can be used, you can install the released version of MaOEA from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("MaOEA")
```

