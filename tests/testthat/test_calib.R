test_that("calib='affine' option in 'vsn2' works with 'predict' method for resulting 'vsn' object.", {
  # motivated by https://support.bioconductor.org/p/92235/
  nr   <- 500
  nc   <- 10
  dat  <- exp(matrix(rnorm(nr*nc), nr, nc) + log(seq_len(nr)))
  fit1 <- vsn2(dat)
  fit2 <- vsn2(dat, calib = "none")
  vs1  <- predict(fit1, newdata = dat)
  vs2  <- predict(fit2, newdata = dat)
})


