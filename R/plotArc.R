#' Plots parabolic arc
#'
#' @param from,to x coordinates of arc
#' @param top highest point of arc
#' @param n number of points
#' @param y.base bottom coordinate of arc
#' @param ... other parameters of lines functoin
plotArc = function(from,to,top,n=100,y.base=0,...){
  len = to - from
  x = seq(from=0,to=len,length.out = n)
  y = x*4*top/len - x^2*(4*top/len^2)
  lines(x+from,y+y.base,...)
}