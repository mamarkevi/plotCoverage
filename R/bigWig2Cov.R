#'
#' @param bw annotation file
#'
#' @export
bigWig2Cov = function(bw){
  bw = as.data.frame(bw)
  bw = bw[order(bw$start),]
  start = bw$start[1]
  stop = bw$end[nrow(bw)]
  r = rep(0,stop-start+1)
  for(i in 1:nrow(bw)){
    r[(bw$start[i]:bw$end[i])-start+1] = bw$score[i]
  }
  list(cov=r,x=start:stop)
}
