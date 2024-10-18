#' @param sid
#' @param jxn
#' @param gene.grange
#'
#' @export
getRecountCov = function(sid,jxn,gene.grange){
  require(recount3)
  require(rtracklayer)

  bigwig = jxn$BigWigURL[colnames(jxn)==sid]

  bw = rtracklayer::import.bw(bigwig,which=gene.grange)
  r = bigWig2Cov(bw)
  j = jxn[,sid]
  jxns = findOverlaps(jxn@rowRanges,gene.grange)
  j = jxn[jxns@from,]
  #introns can be duplicated...
  nms=unique(j@rowRanges@ranges@NAMES)
  j = j[nms,]

  r$juncs =  cbind(as.data.frame(j@rowRanges)[,c('start','end','strand')],score=j@assays@data$counts[,sid])
  r
}
