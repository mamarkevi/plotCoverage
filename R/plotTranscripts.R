plotTranscripts = function(a,
                           ylim=c(0,length(unique(a$transcript_id))),
                           xlim=c(ifelse(a$strand[1]=='+',min(a$start),max(a$stop)),ifelse(a$strand[1]=='+',max(a$stop),min(a$start))),
                           xlab=a$chr_id[1],
                           new=TRUE,yspace=0.8,exon.col='black',cds.col='black',
                           text.cex = 0.7,
                           ...){
  
  if(!is.na(exon.col))
    a$exon.col = exon.col
  if(!is.na(cds.col))
    a$cds.col = cds.col
  
  transc = split(a,a$transcript_id)
  transc = transc[order(sapply(transc,function(x){max(x$stop)-min(x$start)}))]
  if(new)
    plot(1,t='n',xlim=xlim,ylim=ylim,yaxt='n',ylab='',xlab=xlab,...)
  ystep = (ylim[2]-ylim[1])/length(transc)
  
  for(i in 1:length(transc)){
    y = ylim[1] + ystep*i - ystep/2
    t = transc[[i]]
    lines(c(min(t$start),max(t$stop)),c(y,y))
    f = t$feature == 'exon'
    if(sum(f)>0)
      rect(t$start[f],y-ystep/2*yspace,t$stop[f],y+ystep/2*yspace,col = 'white',border = t$exon.col[f])
    f = t$feature == 'CDS'
    if(sum(f)>0)
      rect(t$start[f],y-ystep/2*yspace,t$stop[f],y+ystep/2*yspace,col = t$cds.col[f],border = t$cds.col[f])
  }
  text(par('usr')[1],seq(ylim[1]+ystep/2,by = ystep,length.out = length(transc)),sapply(transc,function(x)x$transcript_name[1]),
       adj = c(1,0.5),xpd=T,cex=text.cex)
}