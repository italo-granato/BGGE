# BLMMD

A `R` package to prepare datasets and fit GxE genomic models

## Installation

```R
install.packages("devtools")
devtools::install_github("italo-granato/BLMMD")
```

## Simple use

```R
library(BGLR)
data(wheat)

K <- list(G = list(K = tcrossprod(X)/ncol(X)))
y <- as.vector(wheat.Y[,1])

fit <- BLMMD(y = y, K = K)
```

### Others params

| params |      |
|--------|------|
| XF      | Design matrix for fixed effects |
| ite     | Number of iterations. |
| burn    | Number of iterations to be discarded as burn-in. |
|thin     | Thinin interval. |
| verbose | Should report be printed on screen |
| me      | tolerance for zero. Default is 1e-10 |

## Authors
Italo Granato
Jaime D. Cuevas D.
Francisco J. Luna-Vázquez
