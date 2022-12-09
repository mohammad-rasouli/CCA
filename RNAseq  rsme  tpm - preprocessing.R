library(edgeR)
library(genefilter)

#Import dataset
df <- read.delim("C:/Users/VistaRayaneh/Desktop/mrkxxz/thesis/aban/CCLE_RNAseq_rsem_transcripts_tpm_20180929.txt")
counts <- df[,-c(1:2)]
counts <- as.matrix(counts)
rownames(counts) <- df$gene_id

#CPM cutoff
selectGenes <- function(counts, min.count=1, N=0.95){
  
  lib.size <- colSums(counts)
  MedianLibSize <- median(lib.size)
  CPM.Cutoff <- min.count / MedianLibSize*1e6
  CPM <- edgeR::cpm(counts,lib.size=lib.size)
  
  min.samples <- round(N * ncol(counts))
  
  f1 <- genefilter::kOverA(min.samples, CPM.Cutoff)
  flist <- genefilter::filterfun(f1)
  keep <- genefilter::genefilter(CPM, flist)
  
  ## the same as:
  #keep <- apply(CPM, 1, function(x, n = min.samples){
  #  t = sum(x >= CPM.Cutoff) >= n
  #  t
  #})
  
  return(keep)
}

keep.exprs <- selectGenes(counts, min.count=1, N=0.95)
myFilt <- counts[keep.exprs,]
dim(myFilt)
#Export filtered data set
write.csv(myFilt,"C:/Users/VistaRayaneh/Desktop/myfil.csv")
myfilt <- read.csv("C:/Users/VistaRayaneh/Desktop/myfil.csv")

#SD threshold by 2 method: 
              #1-rowSDs

mds=matrixStats::rowSds(myfilt) # calculate standard deviation of CpGs
ms7 <- myfilt[mds < 7 ,]
ms6 <- myfilt[mds < 6 ,]

              #2-remove_sd_outlier
library(dataPreparation)
sd6 <- remove_sd_outlier(myfilt, n_sigmas = 6, verbose = TRUE) #n_sigmas <=3 = ob <10
sd7 <- remove_sd_outlier(myfilt, n_sigmas = 7, verbose = TRUE)

#Visualization:
pdf(file="Plots.pdf")
#Just CPM cutoff
hist(log2(as.matrix(myfilt[,-1])+1),ylab="",las=2,main="CPM cutoff")

#CMP + SD:
#First method
hist(log2(as.matrix(ms7[,-1])+1),ylab="",las=2,main="CPM + SD (rowSDs / SD < 7) ", col = "cornflowerblue")
hist(log2(as.matrix(ms6[,-1])+1),ylab="",las=2,main="CPM + SD (rowSDs / SD < 6) ", col = "cornflowerblue")
#Second method
hist(log2(as.matrix(sd6[,-1])+1),ylab="",las=2,main="CPM + SD (n_sigmas = 6)", col = "orange")
hist(log2(as.matrix(sd7[,-1])+1),ylab="",las=2,main="CPM + SD (n_sigmas = 7) ", col = "orange")


dev.off()





