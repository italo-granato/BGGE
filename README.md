# BGGE

A `R` package to prepare datasets and fit GxE genomic models

## Installation
To complete installation of BGGE, you have to install a few packages first.

```R
install.packages("devtools")
devtools::install_github('italo-granato/BGGE')

```

## Simple use

```R
library(BGGE)
library(BGLR)
data(wheat)

K <- list(G = list(K = tcrossprod(X)/ncol(X)))
y <- as.vector(wheat.Y[,1])

fit <- BGLR(y = y, K = K)
```

### Others params

| params |      |
|--------|------|
| XF      | Design matrix for fixed effects |
| ite     | Number of iterations. |
| ne      | Number of subjects by environment. |
| burn    | Number of iterations to be discarded as burn-in. |
|thin     | Thinin interval. |
| verbose | Should report be printed on screen? |
| tol      | tolerance for zero. Default is 1e-10 |

## Authors
Italo Granato
Jaime D. Cuevas D.
Francisco J. Luna-Vázquez
