
<p align="center">

<a href="https://github.com/italo-granato/BGGE">
<img src="Logo.png" alt="Bayesian Genomic Linear Models Applied to GE Genome Selection"/>
</a>

<h4 align="center">

Bayesian Genomic Linear Models Applied to GE Genome Selection -
Development version 0.6.2

</h4>

<h4 align="center">

\[Last README update: 2018-05-16\]

</h4>

<p align="center">

[![Release](http://www.r-pkg.org/badges/version-ago/BGGE
"IBCF.MTME release")](https://cran.r-project.org/package=BGGE "CRAN Page")
[![License: GPL
v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg
"LGPL, Version 2.0")](https://www.gnu.org/licenses/gpl-3.0 "LGPL, Version 2.0")
[![Project Status:
Active](http://www.repostatus.org/badges/latest/active.svg
"status")](http://www.repostatus.org/#active "status - active")
[![Downloads](https://cranlogs.r-pkg.org/badges/BGGE
"IBCF.MTME cranlogs")](https://cran.r-project.org/package=BGGE "CRAN Page")

</p>

</p>

# Table of contents

  - [NEWS](#news)
  - [Instructions](#instructions)
      - [Installation](#install)
      - [Load the package](#package)
      - [Quick use](#example1)
      - [Other Params](#params)
  - [How to cite this package](#cite)
  - [Contributions](#contributions)
  - [Authors](#authors)

<h2 id="news">

News of this version (0.6.2)

</h2>

  - Minor improvements

See the last updates in [NEWS](NEWS.md).

<h2 id="instructions">

Instructions for proper implementation

</h2>

<h3 id="install">

Installation

</h3>

To complete installation of dev version of BGGE from GitHub, you must
have previously installed the devtools package.

``` r
install.packages('devtools')
devtools::install_github('italo/BGGE')
```

If you want to use the stable version of BGGE package, install it from
CRAN.

``` r
install.packages('BGGE')
```

<h3 id="package">

Load the package

</h3>

``` r
library(BGGE)
```

<h3 id="example1">

Example of simple usage of the package

</h3>

``` r
library(BGGE)
library(BGLR)
data(wheat)

X <- scale(wheat.X, scale = TRUE, center = TRUE)
rownames(X) <- 1:599
pheno_geno <- data.frame(env = gl(n = 4, k = 599), 
                         GID = gl(n=599, k=1, length = 599*4),
                         value = as.vector(wheat.Y))

# Creating kernel for GE model

K <- getK(Y = pheno_geno, X = X, kernel = "GB", model = "MM")
y <- as.vector(wheat.Y)

fit <- BGGE(y = y, K = K, ne = rep(599, 4))

plot(fit$yHat, y)
```

![](README_files/figure-gfm/BGGEUse-1.png)<!-- -->

<h3 id="params">

Others params

</h3>

| params  |                                                  |
| ------- | ------------------------------------------------ |
| XF      | Design matrix for fixed effects.                 |
| ite     | Number of iterations.                            |
| ne      | Number of subjects by environment.               |
| burn    | Number of iterations to be discarded as burn-in. |
| thin    | Thinin interval.                                 |
| verbose | Should report be printed on screen?              |
| tol     | tolerance for zero. Default is 1e-10             |

<h2 id="cite">

Citation

</h2>

First option, by the paper.

(Comming soon)

Second option, by the package

``` r
citation('BGGE')
```

    ## 
    ## To cite package 'BGGE' in publications use:
    ## 
    ##   Italo Granato, Luna-Vázquez Francisco J. and Cuevas Jaime
    ##   (2018). BGGE: Bayesian Genomic Linear Models Applied to GE
    ##   Genome Selection. R package version 0.6.2.
    ##   https://CRAN.R-project.org/package=BGGE
    ## 
    ## A BibTeX entry for LaTeX users is
    ## 
    ##   @Manual{,
    ##     title = {BGGE: Bayesian Genomic Linear Models Applied to GE Genome Selection},
    ##     author = {Italo Granato and Luna-Vázquez {Francisco J.} and Cuevas Jaime},
    ##     year = {2018},
    ##     note = {R package version 0.6.2},
    ##     url = {https://CRAN.R-project.org/package=BGGE},
    ##   }

<h2 id="contributions">

Contributions

</h2>

If you have any suggestions or feedback, I would love to hear about it.
Feel free to report new issues in [this
link](https://github.com/italo/BGGE/issues/new), also if you want to
request a feature/report a bug, or make a pull request if you can
contribute.

<h2 id="authors">

Authors

</h2>

  - Italo Granato (Author, Maintainer)
  - Jaime D. Cuevas D. (Author)
  - Francisco J. Luna-Vázquez (Author, Maintainer)
