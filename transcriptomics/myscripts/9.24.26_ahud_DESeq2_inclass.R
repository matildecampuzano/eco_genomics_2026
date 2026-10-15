###This is my .r script for exploring  the A. hudsonica RNAseq data
###using DESeq2

## Set your working directory
setwd("~/Projects/eco_genomics_2026/transcriptomics")

## Import the libraries that we're likely to need in this session

library(DESeq2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(ggpubr)
library(vsn)  
library("pheatmap")
library("vsn")

####################################################

### Import our data

####################################################


# Import the counts matrix
countsTable <- read.table("mydata/salmon.isoform.counts.matrix.filteredAssembly", header=TRUE, row.names=1)
#puts a table of the data into environment so we can open
head(countsTable)
dim(countsTable)
#show dimensions of table

countsTableRound <- round(countsTable) # bc DESeq2 doesn't like decimals (and Salmon outputs data with decimals)
head(countsTableRound)
##rounded the numbers in graph and then show head

#import the sample description table
conds <- read.delim("mydata/ahud_samples_R.txt", header=TRUE, stringsAsFactors = TRUE, row.names=1)
head(conds)


####################################################

### Explore data distributions

####################################################

# Let's see how many reads we have from each sample
colSums(countsTableRound)
mean(colSums(countsTableRound))

barplot(colSums(countsTableRound), names.arg=colnames(countsTableRound),cex.names=0.5, las=3,ylim=c(0,21000000))
abline(h=mean(colSums(countsTableRound)), col="blue", lwd=2)
## shows plot of reads per samples

# the average number of counts per gene
rowSums(countsTableRound)
mean(rowSums(countsTableRound)) # [1] 8217.81
median(rowSums(countsTableRound)) # [1] 377

apply(countsTableRound,2,mean) # 2 in the apply function does the action across columns
apply(countsTableRound,1,mean) # 1 in the apply function does the action across rows
hist(apply(countsTableRound,1,mean),xlim=c(0,10000), ylim=c(0,60000),breaks=10000)

#Define our model and create a DESeq2 object

####################################################

### Start working with DESeq2!

####################################################

#### Create a DESeq object and define the experimental design here with the tilda

dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData=conds, 
                              design= ~ generation + treatment)

dim(dds)

# Filter out genes with too few reads - remove all genes with counts < 15 in more than 75% of samples, so ~28)
## suggested by WGCNA on RNAseq FAQ

dds <- dds[rowSums(counts(dds) >= 15) >= 28,]
nrow(dds) 
# [1] 25260, that have at least 15 reads (a.k.a counts) in 75% of the samples

# Run the DESeq model to test for differential gene expression
dds <- DESeq(dds)

# List the results you've generated

resultsNames(dds)

# Copy the results names:
# [1] "Intercept"            "generation_F11_vs_F0" "generation_F2_vs_F0" 
# [4] "generation_F4_vs_F0"  "treatment_OA_vs_AM"   "treatment_OW_vs_AM"  
# [7] "treatment_OWA_vs_AM" 

## Now lets transform data so we can make a PCA plot!!
# The goal of transformation "is to remove the dependence of the variance 
#on the mean, particularly the high variance of the logarithm of count data 
# when the mean is low."

# this gives log2(n + 1)
ntd <- normTransform(dds)
meanSdPlot(assay(ntd))

# Variance stabilizing transformation
vsd <- vst(dds, blind=FALSE)
meanSdPlot(assay(vsd))

##Look for any outliers in samples

sampleDists <- dist(t(assay(vsd)))

library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$line, vsd$generation, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

# Note any outliers - maybe AH (OA) F2 rep2, prolly because its small amount

sampleTree <- hclust(dist(sampleDists), method="average")
# plot
plot(sampleTree, main="Sample clustering to detect outliers", sub="", xlab="",cex.lab=1.5, cex.axis=1.5, cex.main=2)


# first transform the data for plotting using variance stabilization
vsd <- vst(dds, blind=FALSE)

pcaData <- plotPCA(vsd, intgroup=c("treatment","generation"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData,"percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=treatment, shape=generation)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

###############################################################
## make it prettier
# Let's plot the PCA by generation in four panels to make it PRETTYYYY

data <- plotPCA(vsd, intgroup=c("treatment","generation"), returnData=TRUE)
percentVar <- round(100 * attr(data,"percentVar"))

###########  

dataF0F4F11 <- subset(data, generation %in% c("F0", "F4", "F11"))

PCA_AM_OWA <- ggplot(dataF0F4F11, aes(PC1, PC2)) +
  geom_point(size=7, stroke = 1.5, aes(fill=treatment, shape=generation)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ylim(-25, 30) + xlim(-30, 65) +
  scale_shape_manual(values=c(21,22,24), labels = c("F0", "F4", "F11"))+
  scale_fill_manual(values=c('#6699CC', "#CC3333"), labels = c("Ambient", "OWA"))+
  theme_bw() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 4))+
  theme(text = element_text(size = 20)) +
  theme(legend.title = element_blank())

PCA_AM_OWA
