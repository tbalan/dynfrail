
<!-- README.md is generated from README.Rmd. Please edit that file -->
This is an R package for fitting semiparametric dynamic frailty models with the EM algorithm. The hazard for individual *j* from cluster *i* is specified as:
*λ*<sub>*i**j*</sub>(*t*|*Z*<sub>*i*</sub>(*t*)) = *Z*<sub>*i*</sub>(*t*)exp(*β*<sup>⊤</sup>*x*<sub>*i**j*</sub>(*t*))*λ*<sub>0</sub>(*t*).
 The model used here is described in detail in [Putter & van Houwelingen (2015)](https://doi.org/10.1093/biostatistics/kxv002). The distribution of *Z*<sub>*i*</sub>(*t*) is described by two parameters: *θ*, that is an inverse-variability parameter of *Z*<sub>*i*</sub>(*t*) for a fixed *t*, and *λ*, that describes the autocorrelation of the process, so that for *t*<sub>1</sub> ≤ *t*<sub>2</sub>
*c**o**r*(*Z*<sub>*i*</sub>(*t*<sub>1</sub>),*Z*<sub>*i*</sub>(*t*<sub>2</sub>)) = exp(*λ*(*t*<sub>2</sub> − *t*<sub>1</sub>)).

The estimation process is that for fixed (*θ*, *λ*) the maximized profile likelihood is calculated, i.e. maximized with respect to (*β*, *λ*<sub>0</sub>). This profile likelihood is finally maximized itself.

Installation
------------

The development version from `GitHub`:

``` r
devtools::install_github("tbalan/dynfrail")
```

The following packages are needed to build `dynfrail`:

``` r
install.packages(c("RcppArmadillo", "tibble", "magrittr", "dplyr", "tidyr"))
```

The functioning of the package is described in the documentation of the main fitting function, `dynfrail()`.

Features
--------

-   gamma, PVF, compount Poisson, inverse Gaussian distributions
-   flexible adjustment of estimation parameters
-   semiparametric *Z*(*t*) that changes values at every *t* or piecewise constant *Z*(*t*)
-   clustered survival data & recurrent events (calendar time or gaptime) ar supported

Functions
---------

-   `dynfrail()` has a friendly syntax very similar to the `frailtyEM` package: next to a `formula` and `data` argument, the `distribution` argument is used to specify the distribution parameters and the `control` parameter is used for controling the precision of the estimation.
-   `dynfrail_prep()` and `dynfrail_fit()` are used internally by `dynfrail()` but are made user-available. The first one prepares the input of `dynfrail()` to make it suitable for the actual EM algorithm. The second one performs one EM algorithm for fixed (*θ*, *λ*) to estimate the maximum (*β*, *λ*<sub>0</sub>).

Limitations
-----------

-   slow even for medium sized data sets. It is recommended to start with a small number of piecewise constant intervals and/or a subset of the data
-   no direct standard errors for (*θ*, *λ*).
