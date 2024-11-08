# Instalation and loading of package
devtools::install_github("mamarkevi/plotCoverage")
library('plotCoverage')

# download genome annotation and test samples from ./example directory and load as here

# loading gtf
gtf = loadEnsGTF('./example/a.gtf')

# Samples of 1st condition 
c1 = c('./example/1.bam')
# Samples of 2nd condition
c2 = c('./example/2.bam')

# 1. Plot coverage for first and second condition separately
#1
cov = getReadCoverage(c1,
                      chr = '7',
                      start = 91980000,
                      end = 91993000)

ymax=max(cov$cov,cov$juncs$score)
plotReadCov(cov,ylim=c(-0.4*ymax,ymax),bty='n',min.junc.cov.f = 0.05,xlab='7',ylab='Coverage',
            xlim = c(cov$start,cov$end)) # ylim to leave space for annotation

plotTranscripts(gtf[gtf$gene_name=='AKAP9',],ylim = c(-0.4*ymax,-0.05*ymax),new = F)

#2
cov = getReadCoverage(c2,
                       chr = '7',
                       start = 91980000,
                       end = 91993000)

ymax=max(cov$cov,cov$juncs$score)
plotReadCov(cov,ylim=c(-0.4*ymax,ymax),bty='n',min.junc.cov.f = 0.05,xlab='7',ylab='Coverage',
            xlim = c(cov$start,cov$end)) # ylim to leave space for annotation
# transcripts as 
plotTranscripts(gtf[gtf$gene_name=='AKAP9',],ylim = c(-0.4*ymax,-0.05*ymax),new = F)


# 2. Plot two conditions together - use for comparison of exon coverage
par(mfrow=c(3,1))
#1
cov = getReadCoverage(c1,
                      chr = '7',
                      start = 91980000,
                      end = 91993000)

ymax=max(cov$cov,cov$juncs$score)
plotReadCov(cov,ylim=c(-0.4*ymax,ymax),bty='n',min.junc.cov.f = 0.05,xlab='7',ylab='Coverage',
            xlim = c(cov$start,cov$end)) # ylim to leave space for annotation

#2
cov = getReadCoverage(c2,
                      chr = '7',
                      start = 91980000,
                      end = 91993000)

ymax=max(cov$cov,cov$juncs$score)
plotReadCov(cov,ylim=c(-0.4*ymax,ymax),bty='n',min.junc.cov.f = 0.05,xlab='7',ylab='Coverage',
            xlim = c(cov$start,cov$end)) # ylim to leave space for annotation
# transcripts as single plot 
plotTranscripts(gtf[gtf$gene_name=='AKAP9',],xlim=c(cov$start,cov$end),ylim = c(-0.4*ymax,-0.05*ymax),new = T)

dev.off()







