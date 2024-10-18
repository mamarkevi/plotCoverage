#' @param f file with ensemble annotation
#' @param features
#'
#' @export

loadEnsGTF = function(f,features=NULL){
  r = read.table(f,sep='\t')
  if(!is.null(features ))
    r = r[r$V3 %in% features,]
  a = lapply(strsplit(r$V9,';\\s?',perl=T),function(x){x=strsplit(x,'[ =]');setNames(sapply(x,'[',2),sapply(x,'[',1))})
  names = unique(unlist(lapply(a,names)))
  a = do.call(rbind,lapply(a,'[',names))
  colnames(a) = names
  r = r[,c(1:5,7)]
  colnames(r) = c('chr_id','type','feature','start','stop','strand')
  cbind(r,a)
}
