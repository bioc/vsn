##------------------------------------------------------------
## a wrapper for vsn to be used as a normalization method in
## the package affy
## D P Kreil <bioc07@kreil.org>
##------------------------------------------------------------
normalize.AffyBatch.vsn=function(abatch,
                                 reference,
                                 strata=NULL,
                                 subsample=30000L,
                                 subset,
                                 log2scale=TRUE,
                                 log2asymp=FALSE,
                                 ...) {

  if(is.na(log2scale)||is.na(log2asymp)||(log2scale&&log2asymp))
    stop("'log2asymp' and 'log2scale' must not both be TRUE, and not be NA.")
  
  ind = indexProbes(abatch,"pm")
  if(!missing(subset))
    ind = ind[subset]
  ind = unlist(ind)
       
  vsn2res = vsn2(intensity(abatch)[ind,],reference=reference,
                 returnData=FALSE,subsample=subsample,...)
  
  description(abatch)@preprocessing = c(description(abatch)@preprocessing, 
                       list(vsnReference=vsn2res))
  
  trsfx = predict(vsn2res,newdata=intensity(abatch),log2scale=log2scale)

  ## irrelevant affine transformation for ~ log2(x) for x>>1
  if (log2asymp)
    trsfx=(trsfx-log(2))/log(2)
   
  intensity(abatch)<-2^trsfx
  return(abatch)
}
