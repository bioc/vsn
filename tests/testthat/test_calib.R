test_that("Make sure that if the calib='affine' option is used in 'vsn2', the 'predict' method works on the resulting 'vsn' object.", {
  ## motivated by https://support.bioconductor.org/p/92235/
  set.seed(0xbadbeef)  
  nr   <- 512
  nc   <- 8
  dat  <- exp(matrix(rnorm(nr*nc), nr, nc) + log(seq_len(nr)))
  fit1 <- vsn2(dat)
  fit2 <- vsn2(dat, calib = "none")
  vs1  <- predict(fit1, newdata = dat)
  vs2  <- predict(fit2, newdata = dat)
})


