library("vsn")
data("kidney")

pd =  pData(kidney)
ex =  exprs(kidney)

rownames(pd) = colnames(ex)

vmd = data.frame(labelDescription=I("The scanner channel Cy3 or Cy5"))
rownames(vmd) = colnames(pd)
ad = new("AnnotatedDataFrame", data=pd, varMetadata=vmd)
    
kidney = new("ExpressionSet", exprs=ex, phenoData=ad)

save(kidney, file="../../data/kidney.RData")
