library(affy)
library(vsn)
library(affydata)
data(Dilution)
normalize.AffyBatch.methods <- c(normalize.AffyBatch.methods, "vsn")

print(system.time( {
  es1 = express(Dilution, normalize.method = "vsn",  niter = 4, bg.correct = FALSE)
} ))

print(system.time( {
  es2 = express(Dilution)
} ))
