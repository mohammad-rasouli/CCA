library(DESeq2)
library(ggplot2)


# load counts table 
tbl <- read.csv("Count.csv")

# load gene annotations 
annot <- read.csv('annot.csv')
rownames(annot) <- annot$GeneID


## Selecting sample groups (T47D monoculture and T47D+HS21 co-culture)

sample <- tbl[,c('GSM5870317', 'GSM5870318', 'GSM5870319', "GSM5870347", "GSM5870370",
                 "GSM5870367", "GSM5870368", "GSM5870369", "GSM5870376", "GSM5870378"
                )]

gs <- factor(c(0,0,0,0,0, 1,1,1,1,1))
groups <- make.names(c("T47D mono","T47D-HS21a Co-culture"))
levels(gs) <- groups

sample_info <- data.frame(Group = gs, row.names = colnames(sample))

# Filtering low expressed genes
keep = rowSums(sample >= 10) >= min(table(gs))
tbl <- sample[keep, ]

# Identifying differential expressed genes by DESeq2
ds <- DESeqDataSetFromMatrix(count = tbl, colData = sample_info, design = ~Group)
ds <- DESeq(ds, test="Wald", sfType="poscount")

# Extract results for top genes table
r <- results(ds, contrast=c("Group", groups[1], groups[2]), alpha=0.05, pAdjustMethod ="fdr")

# Filter based on padj and log2FoldChange criteria
filtered_genes <- r[r$padj < 0.05 & abs(r$log2FoldChange) >= 1,]

tT <- merge(as.data.frame(filtered_genes), annot, by=0, sort=F)

tT <- subset(tT, select=c("GeneID","padj","pvalue","lfcSE","stat","log2FoldChange","baseMean","Symbol","Description"))

# Print tT dataframe with filtered results
print(tT)
setwd("C:/Users/mohammad/OneDrive/Desktop/In vitro")
# Save tT as a CSV file
write.csv(tT, file = "T47d vs T from Coculture org.csv", row.names = FALSE)




#Visuallizing

#Volcano plot
library(EnhancedVolcano)


# Create a volcano plot using EnhancedVolcano
EnhancedVolcano(tT,
                lab = tT$Symbol,
                x = "log2FoldChange",
                y = "padj",
                title = "Volcano Plot of Differentially Expressed Genes",
                pCutoff = 0.05,
                FCcutoff = 1,
                xlim = c(-10, 10))


#*********************Control Groups*****************************************************
## T47D and HS27a from T47D+HS27a Co-Culture 

sample <- tbl[,c('GSM5870317', 'GSM5870318', 'GSM5870319', "GSM5870347", "GSM5870370",
                  "GSM5870364", "GSM5870365", "GSM5870366", "GSM5870375", "GSM5870377")]


gs <- factor(c( 0,0,0,0,0,1,1,1,1,1))
groups <- make.names(c("T47D", "H from T+H"))
levels(gs) <- groups

sample_info <- data.frame(Group = gs, row.names = colnames(sample))

# Filtering low expressed genes:
keep = rowSums(sample >= 10) >= min(table(gs))
tbl <- sample[keep, ]

# DESeq2
ds <- DESeqDataSetFromMatrix(count = tbl, colData = sample_info, design = ~Group)
ds <- DESeq(ds, test="Wald", sfType="poscount")

# Extract results for top genes table
r <- results(ds, contrast=c("Group", groups[1], groups[2]), alpha=0.05, pAdjustMethod ="fdr")

# Filter based on padj and log2FoldChange criteria
filtered_genes <- r[!is.na(r$padj) & !is.na(r$log2FoldChange) & r$padj < 0.05 & abs(r$log2FoldChange) >=2,]


#filtered_genes <- r[r$padj < 0.05 & abs(r$log2FoldChange) >= 1,]



#tT <- r[order(-abs(r$log2FoldChange))[1:300],]
tT <- merge(as.data.frame(filtered_genes), annot, by=0, sort=F)


tT <- subset(tT, select=c("EnsemblGeneID", "GeneID","padj","pvalue","lfcSE","stat","log2FoldChange","baseMean"))

# Print tT dataframe with filtered results
print(tT)
setwd("C:/Users/mohammad/OneDrive/Desktop/In vitro")

# Save tT as a CSV file
write.csv(tT, file = "T47d vsH from Coculture org.csv", row.names = FALSE)


## T47D vs HS27A mono culture

sample <- tbl[,c('GSM5870317', 'GSM5870318', 'GSM5870319', "GSM5870347", "GSM5870370",
                 "GSM5870309", "GSM5870323", "GSM5870324", "GSM5870325", "GSM5870372")]


gs <- factor(c( 0,0,0,0,0,1,1,1,1,1))
groups <- make.names(c("T47D", "HS27A"))
levels(gs) <- groups

sample_info <- data.frame(Group = gs, row.names = colnames(sample))

# Filtering low expressed genes:
keep = rowSums(sample >= 10) >= min(table(gs))
tbl <- sample[keep, ]

# DESeq2
ds <- DESeqDataSetFromMatrix(count = tbl, colData = sample_info, design = ~Group)
ds <- DESeq(ds, test="Wald", sfType="poscount")

# Extract results for top genes table
r <- results(ds, contrast=c("Group", groups[1], groups[2]), alpha=0.05, pAdjustMethod ="fdr")

# Filter based on padj and log2FoldChange criteria
filtered_genes <- r[!is.na(r$padj) & !is.na(r$log2FoldChange) & r$padj < 0.05 & abs(r$log2FoldChange) >=2,]


#filtered_genes <- r[r$padj < 0.05 & abs(r$log2FoldChange) >= 1,]



#tT <- r[order(-abs(r$log2FoldChange))[1:300],]
tT <- merge(as.data.frame(filtered_genes), annot, by=0, sort=F)


tT <- subset(tT, select=c("EnsemblGeneID", "GeneID","padj","pvalue","lfcSE","stat","log2FoldChange","baseMean"))

# Print tT dataframe with filtered results
print(tT)
setwd("C:/Users/mohammad/OneDrive/Desktop/In vitro")

# Save tT as a CSV file
write.csv(tT, file = "T47d vs HS27a mono.csv", row.names = FALSE)

