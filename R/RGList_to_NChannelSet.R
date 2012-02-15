## I copied this file from https://hedgehog.fhcrc.org/bioconductor/trunk/madman/Rpacks/convert/R/RGList_to_NChannelSet.R
## on 15.2.2012. The subversion long stated for the file of origin:
##     "Last Changed Rev: 36027, Last Changed Date: 2008-12-16 22:48:47 +0100 (Tue, 16 Dec 2008)"


setAs("RGList", "NChannelSet", function(from) {

  ## assayData
  assayData <- with(from, {
    if (!exists("other", inherits=FALSE)) {
      elts <- list(R=R, G=G)
    } else {
      if (is.null(names(other)) ||
          !all(sapply(names(other), nzchar)))
        stop(paste("RGList 'other' elements must be named, found '",
                   paste(names(other), collapse="', '"),
                   "'", sep=""))
      bad <- names(other) %in% c("R", "G", "Rb", "Gb")
      if (any(bad))
        stop(paste("RGList 'other' elements contain reserved names '",
                   paste(names(other)[bad],
                         collapse="', '"),
                   "'", sep=""))
      
      elts <- c(R=R, G=G, other)
    }
    if (exists("Rb", inherits=FALSE))
      elts[["Rb"]] <- Rb
    if (exists("Gb", inherits=FALSE))
      elts[["Gb"]] <- Gb
    do.call(assayDataNew,
            c(storage.mode="lockedEnvironment", elts))
  })
  
  ## phenoData
  phenoData <- 
    if (!is.null(from$targets)) {
      new("AnnotatedDataFrame", data=from$targets)
    } else {
      new("AnnotatedDataFrame",
          data=data.frame(rep(0, ncol(from)))[,FALSE])
    }
  ## featureData
  fData <-
    if (!is.null(from$genes))
      from$genes
    else
      data.frame(x=rep(0,nrow(from)))[,FALSE]
  if (!is.null(rownames(assayData[["R"]])))
    row.names(fData) <- rownames(assayData[["R"]])
  featureData <- new("AnnotatedDataFrame", data = fData)
  if (!is.null(from$weights))
    if ("weights" %in% names(df))
      warning("RGList 'genes' contains column 'weights'; 'wt.fun' weights discarded")
    else
      featureData[["weights",
                   labelDescription="calculated, from RGList"]] <-
                     from$weights
  ## experimentData
  other <- 
    if (!is.null(from$source))
      list(source=from$source)
    else
      list()
  experimentData <- new("MIAME",
                        other=c("converted from marrayRaw", other))
  new("NChannelSet",
      assayData = assayData,
      featureData = featureData,
      phenoData = phenoData,
      experimentData=experimentData)
})
