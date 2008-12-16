## Wolfgang Huber 2007-2008
## 
## Based on code by Martin Morgan (27 August 2007) from the package "convert"

setAs("RGList", "NChannelSet", function(from) {

  ## assayData
  assayData <- with(from, {
    if (exists("other", inherits=FALSE))
      note("Ignoring slot 'other'.")

    elts <- list(R=R, G=G)

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
