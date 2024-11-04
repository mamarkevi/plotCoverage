#' Plots read coverage
#'
#' @param r read coverage; output of \code{\link{getReadCoverage}}
#' @param min.junc.cov numeric, plots only junctions (introns) with coverage not less than \code{min.junc.cov}
#' @param min.junc.cov.f numeric, plots only junctions (introns) with coverage not less than \code{min.junc.cov.f} of maximal coverage in the region
#' @param plot.junc.only.within logical, plot only junction with both ends within the region, FALSE plots all junctions with at least one end within region. NA plot all junctions overlapping the region.
#' @param ylim,xlim see \code{\link{plot}}
#' @param reverse reverse x coordinates
#' @param junc.col color for junction line. Individual color could be specified for each junction
#' @param junc.lwd line width for junction line
#' @param ... other parameters for plot function
#'
#' @export
plotReadCov = function(r,min.junc.cov=0,min.junc.cov.f=0,plot.junc.only.within=FALSE,ylim=NULL,
                       xlim=range(r$x),reverse=FALSE,junc.col='blue',junc.lwd=3,bottom.mar=0,...){
  f = r$x >= xlim[1] & r$x <=xlim[2]
  r$x = r$x[f]
  r$cov = r$cov[f]
  if(nrow(r$juncs)>0)
    r$juncs$col = junc.col
  r$juncs = r$juncs[r$juncs$start <= xlim[2] & r$juncs$end >=xlim[1] & r$juncs$score >= min.junc.cov,]
  if(!is.na(plot.junc.only.within)){
    if(plot.junc.only.within){
      r$juncs = r$juncs[r$juncs$start > xlim[1] & r$juncs$end < xlim[2],]
    }else{
      r$juncs = r$juncs[(r$juncs$start > xlim[1] & r$juncs$start < xlim[2]) | (r$juncs$end > xlim[1] & r$juncs$end < xlim[2]),]
    }
  }
  
  start = r$x[1]
  end = r$x[length(r$x)]
  r$cov[c(1,length(r$cov))] = 0
  if(is.null(ylim)){
    ylim = c(0,max(r$cov,ifelse(nrow(r$juncs)>0,max(r$juncs$score),1)))
    ylim = c(-bottom.mar*ylim[2],ylim[2])
  }
  r$juncs = r$juncs[r$juncs$score >= min.junc.cov.f * ylim[2],]
  if(reverse)
    xlim=rev(xlim)
  plot(r$x,r$cov,t='n',ylim=ylim,xlim=xlim,yaxt='n',...)
  axis(2,at=c(0,ylim[2]),labels = c('',ylim[2]))
  polygon(r$x,r$cov,col = 'gray',border=NA)
  if(nrow(r$juncs)>0)
    for(i in 1:nrow(r$juncs)){
      start = r$juncs$start[i]
      stop = r$juncs$end[i]
      # to make junction max height always within region
      if(start<xlim[1] & stop >= xlim[2]){
        start = xlim[1] - mean(xlim)
        stop = xlim[2] + mean(xlim)
      }else if (start < xlim[1]){
        start = max(start,xlim[1] - (stop - xlim[1]))
      }else if (stop > xlim[2]){
        stop = min(stop,xlim[2] + (xlim[2]-start))
      }
      plotArc(start,stop,r$juncs$score[i],col=r$juncs$col[i],lwd=junc.lwd)
    }
  invisible(ylim)
}
