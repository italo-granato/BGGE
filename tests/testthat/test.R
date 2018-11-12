library(testthat)
library(BGGE)
library(BGLR)

context('Kernel - With One Environment')
test_that('Generate Kernel with BGLR data', {
  data(wheat)
  
  X <- scale(wheat.X, scale = TRUE, center = TRUE)
  rownames(X) <- 1:599
  pheno_geno <- data.frame(env = 'One', 
                           GID = rownames(X),
                           value = wheat.Y[,1])
  
  # Expect error with multiple environments.
  expect_error(getK(Y = pheno_geno, X = X, kernel = "GB", model = "MM"))
  expect_error(getK(Y = pheno_geno, X = X, kernel = "GK", model = "MM"))
  
  # Test GB Kernel
  GB_Kernel <- getK(Y = pheno_geno, X = X, kernel = "GB", model = "SM")
  expect_output(str(GB_Kernel), 'List of 2')
  expect_equal(names(GB_Kernel), 'G')
  expect_is(GB_Kernel$G$Kernel, 'matrix')
  expect_equal(GB_Kernel$G$Type, 'D')
  
  # Test GK Kernel
  GK_Kernel <- getK(Y = pheno_geno, X = X, kernel = "GK", model = "SM")
  expect_output(str(GK_Kernel), 'List of 2')
  expect_equal(names(GK_Kernel), 'G')
  expect_is(GK_Kernel$G$Kernel, 'matrix')
  expect_equal(GK_Kernel$G$Type, 'D')
})


context('Kernel - Multi Environment')
test_that('Generate Kernel with BGLR data', {
  data(wheat)
  
  X <- scale(wheat.X, scale = TRUE, center = TRUE)
  rownames(X) <- 1:599
  pheno_geno <- data.frame(env = gl(n = 4, k = 599), 
                           GID = gl(n = 599, k = 1, length = 599 * 4),
                           value = as.vector(wheat.Y))
  
  # Expect error with single environments.
  expect_error(getK(Y = pheno_geno, X = X, kernel = "GB", model = "SM"))
  expect_error(getK(Y = pheno_geno, X = X, kernel = "GK", model = "SM"))
  
  # Test GB Kernel
  GB_Kernel <- getK(Y = pheno_geno, X = X, kernel = "GB", model = "MM")
  expect_output(str(GB_Kernel), 'List of 2')
  expect_equal(names(GB_Kernel), 'G')
  expect_is(GB_Kernel$G$Kernel, 'matrix')
  expect_equal(GB_Kernel$G$Type, 'D')
  
  # Test GK Kernel
  GK_Kernel <- getK(Y = pheno_geno, X = X, kernel = "GK", model = "MM")
  expect_output(str(GK_Kernel), 'List of 2')
  expect_equal(names(GK_Kernel), 'G')
  expect_is(GK_Kernel$G$Kernel, 'matrix')
  expect_equal(GK_Kernel$G$Type, 'D')
})

context('Fitting')

test_that('Fit and predict BGLR data', {
  data(wheat)
  
  X <- scale(wheat.X, scale = TRUE, center = TRUE)
  rownames(X) <- 1:599
  pheno_geno <- data.frame(env = gl(n = 4, k = 599), 
                           GID = gl(n = 599, k = 1, length = 599 * 4),
                           value = as.vector(wheat.Y))
  
  # Creating kernel for GE model
  K <- getK(Y = pheno_geno, X = X, kernel = "GB", model = "MM")
  
  fm <- BGGE(y = pheno_geno$value, K = K, ne = rep(599, 4))

  # Test fitted model
  expect_is(fm, 'BGGE')
  expect_output(str(fm), 'List of 9')
  expect_is(fm$y, 'numeric')
  expect_is(fm$yHat, 'numeric')
})
