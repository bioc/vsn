##------------------------------------------------------------
## a wrapper for vsn to be used as a normalization method in
## the package affy
##------------------------------------------------------------
## D P Kreil <bioc07@kreil.org>
normalize.AffyBatch.vsn<-function(abatch,
                                  reference,
                                  strata=NULL,
                                  subsample=30000L,
                                  subset,
                                  log2scale=TRUE,
                                  log2asymp=FALSE,
                                  ...) {
  require(affy);
  if (missing(subset)) {
    ind<-unlist(indexProbes(abatch,"pm"));
  } else {
    ind<-unlist(indexProbes(abatch,"pm")[subset]);
  }
  if (!is.na(subsample)) {
    if (!is.numeric(subsample))
      stop("'subsample' must be numeric.");
    if (length(ind)>subsample)
      ind<-sample(ind,subsample);
  }
  if (length(strata)>0) {  # does this make sense for affy? was at probe level!
    stop("normalize.AffyBatch.vsn: strata not supported");
#    strind<-rep(as.integer(1), length(ind));
#    strata=factor(integer(0), levels="all");
  }
  vsn2res<-vsn2(intensity(abatch)[ind,],reference=reference,
                returnData=FALSE,...);
  description(abatch)@preprocessing = c(description(abatch)@preprocessing, 
                       list(vsnReference=vsn2res));
  trsfx = predict(vsn2res,newdata=intensity(abatch),log2scale=log2scale);
  if (log2asymp) { ## irrelevant affine transformation for ~ log2(x) for x>>1
    if (log2scale) stop("log2asymp==TRUE implies log2scale==FALSE.");
    trsfx<-(trsfx-log(2))/log(2);
  }
  intensity(abatch)<-2^trsfx;
  return(abatch);
}
