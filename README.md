leapgp
================

[![License: GPL
v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![](https://img.shields.io/badge/devel%20version-1.0.0-purple.svg)](https://github.com/knrumsey/leapgp)

<div class="figure">

<img src="inst/logos/LEAPGP.png" alt="This logo was designed by Imagine AI Art Studio" width="40%" />
<p class="caption">
This logo was designed by Imagine AI Art Studio
</p>

</div>

### Description

<!-- README.md is generated from README.Rmd. Please edit that file -->

`leapgp` is an R package designed for fast sequential emulation of
computer models. Sequential emulation is common in applications
requiring MCMC. The `leapgp` package is based on [Rumsey et
al. (2023)](https://onlinelibrary.wiley.com/doi/pdf/10.1002/sta4.576)
and can be viewed as a global model extension of the `laGP` method of
Apley and Gramacy 2015.

To install the package, use

``` r
# install.packages("devtools")
devtools::install_github("knrumsey/leapgp")
```

The method is described in [Rumsey et
al. (2023)](https://onlinelibrary.wiley.com/doi/pdf/10.1002/sta4.576)
and explicit code examples are given in the vignette associated with
this package.

## References

Rumsey, Kellin N., Gabriel Huerta, and J. Derek Tucker. “A localized
ensemble of approximate Gaussian processes for fast sequential
emulation.” Stat (2023): e576.

Gramacy, Robert B., and Daniel W. Apley. “Local Gaussian process
approximation for large computer experiments.” Journal of Computational
and Graphical Statistics 24.2 (2015): 561-578.
