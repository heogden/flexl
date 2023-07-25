
<!-- README.md is generated from README.Rmd. Please edit that file -->

# flexl

<!-- badges: start -->

<!-- badges: end -->

The goal of flexl is to …

## Installation

To improve speed of C++ code, we need to switch off debugging. In
development, can do this by running

``` r
pkgbuild::compile_dll(debug = FALSE)
devtools::load_all()
```

(installing the package with `install.packages` switches off debugging
automatically).

When compiling, I set the Makevars (`~/.R/Makevars`) to

    CPPFLAGS = -w

to switch off the warnings. This was necessary because the code was very
slow to compile if all warnings were printed (because of things included
from Boost and Eigen).

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(flexl)
## basic example code
```
