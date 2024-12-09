#' Returns the summary coverage for region
#' @param l list of coverages for same region to be summed
#'
#' @export

sumCovs = function(l){
  r = l[[1]]
  if(length(l)==1)
    return(r);
  juncs = unique(do.call(rbind,unname(lapply(l,function(c)c$juncs[,1:3]))))
  juncs$score = 0
  juncs[rownames(r$juncs),'score'] = r$juncs$score

  for(i in 2:length(l)){
    if(!all(r$x==l[[i]]$x))
      stop("objects should cover identicall intervals")
    r$cov = r$cov + l[[i]]$cov
    juncs[rownames(l[[i]]$juncs),'score'] = juncs[rownames(l[[i]]$juncs),'score'] + l[[i]]$juncs$score
  }
  r$juncs = juncs
  r
}
