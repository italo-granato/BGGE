# BGGE

A `R` package to prepare datasets and fit GE genomic models

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

X <- scale(wheat.X, scale = TRUE, center = TRUE)
rownames(X) <- 1:599
pheno_geno <- data.frame(env = gl(n = 4, k = 599), 
                         GID = gl(n=599, k=1, length = 599*4),
                         value = as.vector(wheat.Y))

# Creating kernel for GE model

K <- getK(Y = pheno_geno, X = X, kernel = "GB", model = "MM")
y <- as.vector(wheat.Y)

fit <- BGGE(y = y, K = K, ne = rep(599, 4))
```

### Others params

| params |      |
|--------|------|
| XF      | Design matrix for fixed effects. |
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


## Link to dropbox
https://www.dropbox.com/sh/1aezkhn359b3yx8/AAC8XbERcAW5ZPDfHlpspJs0a?dl=0
