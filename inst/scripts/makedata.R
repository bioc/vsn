## ----------------------------------------------------------------------
## This file contains the documentation and R commands that were used to
## prepare the example data in the data subdirectory of the package
## ----------------------------------------------------------------------
dataout   = "../../data"
library(Biobase)

## ------------------------------------------------------------
## 1. lymphoma
## ------------------------------------------------------------
samples   = read.delim("lymphomasamples.txt", as.is=T)
datain    = "/home/whuber/h/VSN/alizadeh"

## CH1  = Cy3 = green =  reference
## CH2  = Cy5 = red   =  sample of interest

nrspots   = 9216
nrsamples = nrow(samples)
qua       = matrix(NA, nrow=nrspots, ncol=2*nrsamples)
pd        = data.frame(name    = I(character(2*nrsamples)),
                       sample  = I(character(2*nrsamples)))
  
for (i in 1:nrsamples) {
  filename = paste(samples$name[i], 'rex.DAT', sep='')
  dat = read.delim(file.path(datain, filename))
  qua[,2*i-1] = dat$CH1I - dat$CH1B
  qua[,2*i]   = dat$CH2I - dat$CH2B
  pd$name[(2*i-1):(2*i)]    = samples$name[i]
  pd$sample[2*i-1]  = "reference"
  pd$sample[2*i]    = samples$sampleid[i]
}
colnames(qua) = pd$sample

lymphoma = new("exprSet",
    exprs = qua,
    phenoData = new("phenoData",
      pData     = pd,
      varLabels = list(name="Name of the Chip", sample="Sample")))
    
save(lymphoma, file=file.path(dataout, 'lymphoma.RData'), compress=TRUE)

## ------------------------------------------------------------
## 2. kidney
## ------------------------------------------------------------
datain = "/net/herkules/raid4/home/whuber/Kidney2"
thehyb = 90
load(file.path(datain, "squa.Rdata"))

dat = (squa[, c("fg.green", "fg.red"), thehyb]
      -squa[, c("bg.green", "bg.red"), thehyb])
rownames(dat) = NULL
colnames(dat) = c("green", "red")

kidney = new("exprSet",
  exprs = dat,
  phenoData = new("phenoData",
    pData = data.frame(channel = c("green", "red")),
    varLabels = list(channel="green: 532 nm, dye=Cy3; red: 635 nm, dye=Cy5")))

save(kidney, file=file.path(dataout, 'kidney.RData'), compress=TRUE)
