## This test makes sure that the file RGList_to_NChannelSet.R in this package
##   is the same as the one in the 'convert' package
## It makes the frivolous assumption that "madman/Rpacks" is in the current
##  path and that the source directories for these packages are
##  found underneath. Let's see for many build systems that is true.

p = strsplit(getwd(), .Platform$file.sep)[[1L]]
w1 = which(p=="madman")
w2 = which(p=="Rpacks")
stopifnot(length(w1)==1L, length(w2)==1L, w2==w1+1L)
rt = paste(p[1L:w2], collapse= .Platform$file.sep)
zup = readLines(file.path(rt, "vsn",     "R", "RGList_to_NChannelSet.R"))
zap = readLines(file.path(rt, "convert", "R", "RGList_to_NChannelSet.R"))
cat(if(identical(zup, zap)) "Bingo" else "Zapperlot", "\n")
