#' Extract read coverage from bam files
#'
#' @param bams character vector with paths to bam-files (bai-files requires to be storage in the same directory)
#' @param chr contig name
#' @param start,end coordinates of region
#' @param strand strand, NA for unstranded (default)
#' @param scanBamFlags list of flags to filter reads, see GenomicAlignments::scanBamFlag
#' @param tagFilter list of tags value to filter by. Can be used to filter cells in single cell bams, for example tagFilter=list(CB=barcodes)
#'
#' @return list with three items: x (chr coordinates); cov - number of reads mapped to chr position; juncs - data.frame with introns
#' @export
getReadCoverage = function(bams,chr,start,end,strand=NA,scanBamFlags=list(),tagFilter=list()){
  require(GenomicAlignments)
  if(start>end){
    t = start
    start = end
    end=t
  }
  r = list(x = start:end,
           cov = NULL,
           juncs=NULL,
           chr=chr,
           start=start,
           end=end)
  scanBamFlags$isMinusStrand = strand==-1
  flags = do.call(scanBamFlag,scanBamFlags)
  param = ScanBamParam(flag=flags,which=GRanges(chr, IRanges(start, end)),tagFilter = tagFilter)
  i = 1
  for(b in bams){
    cat('\r',i,'     ')
    i = i + 1
    bam = readGAlignments(b,param = param)
    cov=coverage(bam)[[chr]][start:end]
    juncs = as.data.frame(summarizeJunctions(bam))
    rownames(juncs)=paste(juncs$seqnames,juncs$start,juncs$end,sep='-')
    if(is.null(r$cov)){
      r$cov=cov
      r$juncs=juncs
    }else{
      r$cov = r$cov + cov
      cmn = intersect(rownames(juncs),rownames(r$juncs))
      r$juncs[cmn,'score'] = r$juncs[cmn,'score'] + juncs[cmn,'score']
      r$juncs = rbind(r$juncs,juncs[setdiff(rownames(juncs),rownames(r$juncs)),])
    }
  }
  invisible(r)
}
