metabolomics <- read.csv("C:/Users/Mohammad/Desktop/new/metabolome.csv")
methylome <- read.csv ("C:/Users/Mohammad/Desktop/new/methyome2.csv")

x <- metabolomics [, 3:227]
y <- methylome [, 3:36]


cc <- cc(x,y)

#import required libraries
library(CCP)
library(CCA)
library(GGally)
library(ggplot2)
##checking the between and within set associations
cormat <- matcor (x,y)
#correlation within x dataset
xcor <- cormat$Xcor
xcor <- as.data.frame(xcor)

#correlation within x dataset
ycor <- cormat$Ycor
ycor <- as.data.frame(ycor)

#correlation between x and y datasets
xycor <- cormat$XYcor
xycor <- as.data.frame(xycor)

#CCA package
cc1 <- cc(x,y)
cc1$cor


#loadings
cc2 <- comput(x, y, cc1)
cc2[3:6]
cc2

# tests of canonical dimensions
rho <- cc1$cor
## Define number of observations, number of variables in first set, and number of variables in the second set.
n <- dim(x)[1]
p <- length(x)
q <- length(y)

## Calculate p-valuesus:
p.asym(rho, n, p, q, tstat = "Wilks")
p.asym(rho, n, p, q, tstat = "Hotelling")
p.asym(rho, n, p, q, tstat = "Pillai")
p.asym(rho, n, p, q, tstat = "Roy")

# standardized x canonical coefficients
s1 <- diag(sqrt(diag(cov(x))))

#for significant CVs
xstd.coeff <- s1 %*% cc1$xcoef[,1:6]
View(as.data.frame(xstd.coeff))

#Sorted Correlation coefficients for X and Y sets

library(dplyr)
library(tidyr)

#For X set(metabolome)

corx1<- 
  cor(x) %>%
  as.data.frame() %>%
  mutate(var1 = rownames(.)) %>%
  gather(var2, value, -var1) %>%
  arrange(desc(value)) %>%
  group_by(value) %>%
  filter(row_number()==1)

#For y set (Methylation)

cory1<- 
  cor(y) %>%
  as.data.frame() %>%
  mutate(var1 = rownames(.)) %>%
  gather(var2, value, -var1) %>%
  arrange(desc(value)) %>%
  group_by(value) %>%
  filter(row_number()==1)

#exporting loadings. as a csv file
write.csv(cc2$corr.X.xscores,"C:/Users/Mohammad/Desktop/X.xscore.csv", row.names = TRUE)
write.csv(cc2$corr.Y.yscores,"C:/Users/Mohammad/Desktop/Y.yscore.csv", row.names = TRUE)

cell_id <- read.csv("C:/Users/Mohammad/Desktop/cell_id.csv")
#creating a separate dataset for methylation CVs
methylCV <- data.frame(row.names = cell_id[,1])

#CV1
methylCV$cv1P <- rowMeans(y[ , c("H3K9me3K14ac0","H3K4me1",
                                 "H3K27me2K36me2","H3K27me1K36me3",
                                 "H3K27me2K36me1")], na.rm=TRUE)

methylCV$cv1N <- rowMeans(y[ , c("H3K9me0K14ac1","H3K9me1K14ac1")], na.rm=TRUE)

#CV2

methylCV$cv2P <- rowMeans(y[ , c("H3K27me1K36me0","H3K9me3K14ac0",
                                 "H3K27me0K36me0","H3K27me0K36me3")], na.rm=TRUE)

methylCV$cv2N <- rowMeans(y[ , c("H3K9me2K14ac1","HH3K9me1K14ac1",
                                 "H3K27me1K36me2","H3K4me1","H3K27ac1K36me1")], na.rm=TRUE)
#CV3

methylCV$cv3P <- rowMeans(y[ , c("H3K9me0K14ac0","H3K27me0K36me0",
                                 "H3K27me1K36me0","H3.3K27me0K36me0")], na.rm=TRUE)

methylCV$cv3N <- rowMeans(y[ , c("H3K27me2K36me2","H3K9me3K14ac1",
                                 "H3K79me2")], na.rm=TRUE)

#CV4
methylCV$cv4P <- rowMeans(y[ , c("H3K27me3K36me1","H3K4me1",
                                 "H3K9me0K14ac0")], na.rm=TRUE)

methylCV$cv4N <- (y[ , "H3K27me1K36me0"])

#CV5
methylCV$cv5P <- (y[ , "H3K9me1K14ac0"])

methylCV$cv5N <- rowMeans(y[ , c("H3K27me1K36me2","H3K27ac1K36me2")], na.rm=TRUE)

#CV6
methylCV$cv6P <- (y[ , ("H3K27ac1K36me0")])

methylCV$cv6N <- (y[ , ("H3K27me0K36me2")])

#CV7
methylCV$cv7P <- rowMeans(y[ , c("H3K9me1K14ac1","H3K27me3K36me1","H3K9me0K14ac1",
                                 "H3K27me3K36me0","H3K27me2K36me1")], na.rm=TRUE)

methylCV$cv7N <- rowMeans(y[ , c("H3K27me0K36me1","H3K27me0K36me2","H3.3K27me0K36me0",
                                 "H3K27me1K36me1")], na.rm=TRUE)


#export 
write.csv(methylCV,"C:/Users/Mohammad/Desktop/methylCV.csv", row.names = TRUE)














