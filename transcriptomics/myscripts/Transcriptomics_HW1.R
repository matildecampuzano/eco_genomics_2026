## Set working directory
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
library(UpSetR)
library(topGO)

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


### Explore data distributions - needed for paper, not code####


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

### Start working with DESeq2!


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


## Make volcano plots for F0,F4,F11####

# Filter 
dds <- dds[rowSums(counts(dds) >= 15) >= 28,]
nrow(dds) 

# Subset the DESeqDataSet to the specific level of the "generation" factor
### Subset to F0, F4, and F11
### Keep only AM and OWA treatments

# F0
dds_F0 <- dds[, dds$generation == "F0" &
                dds$treatment %in% c("AM", "OWA")]

dds_F0$treatment <- droplevels(dds_F0$treatment)
design(dds_F0) <- ~ treatment

# F4
dds_F4 <- dds[, dds$generation == "F4" &
                dds$treatment %in% c("AM", "OWA")]

dds_F4$treatment <- droplevels(dds_F4$treatment)
design(dds_F4) <- ~ treatment

# F11
dds_F11 <- dds[, dds$generation == "F11" &
                 dds$treatment %in% c("AM", "OWA")]

dds_F11$treatment <- droplevels(dds_F11$treatment)
design(dds_F11) <- ~ treatment


### Check the number of samples in each group
table(dds_F0$treatment)
table(dds_F4$treatment)
table(dds_F11$treatment)



# Perform DESeq2 analysis on the subset
dds_F0 <- DESeq(dds_F0)
dds_F4 <- DESeq(dds_F4)
dds_F11 <- DESeq(dds_F11)



### Check on the DE results from the DESeq



resultsNames(dds_F0)
# [1] "Intercept"         "treatment_OWA_vs_AM"

res_OWAvsAM_F0 <- results(dds_F0, name="treatment_OWA_vs_AM", alpha=0.05)
res_OWAvsAM_F0 <- res_OWAvsAM_F0[order(res_OWAvsAM_F0$padj),]
head(res_OWAvsAM_F0)  
summary(res_OWAvsAM_F0)

res_OWAvsAM_F4 <- results(dds_F4, name="treatment_OWA_vs_AM", alpha=0.05)
res_OWAvsAM_F4 <- res_OWAvsAM_F4[order(res_OWAvsAM_F4$padj),]
head(res_OWAvsAM_F4)  
summary(res_OWAvsAM_F4)

res_OWAvsAM_F11 <- results(dds_F11, name="treatment_OWA_vs_AM", alpha=0.05)
res_OWAvsAM_F11 <- res_OWAvsAM_F11[order(res_OWAvsAM_F11$padj),]
head(res_OWAvsAM_F11)  
summary(res_OWAvsAM_F11)

##make volcano plot

volcano_df_F0 <- as.data.frame(res_OWAvsAM_F0)

volcano_df_F0 <- volcano_df_F0 %>%
  mutate(
    sig = case_when(
      padj < 0.05 & log2FoldChange > 1  ~ "Up",
      padj < 0.05 & log2FoldChange < -1 ~ "Down",
      TRUE ~ "NS"
    )
  )

volcano_df_F0 <- ggplot(volcano_df_F0,
                        aes(x = log2FoldChange,
                            y = -log10(padj),
                            color = sig)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c(
    "Down" = "navy",
    "NS"   = "gray",
    "Up"   = "darkred"
  )) +
  geom_vline(xintercept = c(-1, 1),
             linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") +
  theme_classic(base_size = 14) +
  labs(
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value",
    color = NULL,
    title = "OWA vs AM (F0)"
  ) +
  theme(plot.title = element_text(
    hjust = 0.5,
    size = 18,
    face = "bold"
  ))
  
volcano_df_F0

volcano_df_F4 <- as.data.frame(res_OWAvsAM_F4)

volcano_df_F4 <- volcano_df_F4 %>%
  mutate(
    sig = case_when(
      padj < 0.05 & log2FoldChange > 1  ~ "Up",
      padj < 0.05 & log2FoldChange < -1 ~ "Down",
      TRUE ~ "NS"
    )
  )

volcano_df_F4 <- ggplot(volcano_df_F4,
                        aes(x = log2FoldChange,
                            y = -log10(padj),
                            color = sig)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c(
    "Down" = "navy",
    "NS"   = "gray",
    "Up"   = "darkred"
  )) +
  geom_vline(xintercept = c(-1, 1),
             linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") +
  theme_classic(base_size = 14) +
  labs(
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value",
    color = NULL,
    title = "OWA vs AM (F4)"
  ) +
  theme(plot.title = element_text(
    hjust = 0.5,
    size = 18,
    face = "bold"
  ))

volcano_df_F4

volcano_df_F11 <- as.data.frame(res_OWAvsAM_F11)

volcano_df_F11 <- volcano_df_F11 %>%
  mutate(
    sig = case_when(
      padj < 0.05 & log2FoldChange > 1  ~ "Up",
      padj < 0.05 & log2FoldChange < -1 ~ "Down",
      TRUE ~ "NS"
    )
  )

volcano_df_F11 <- ggplot(volcano_df_F11,
                         aes(x = log2FoldChange,
                             y = -log10(padj),
                             color = sig)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c(
    "Down" = "navy",
    "NS"   = "gray",
    "Up"   = "darkred"
  )) +
  geom_vline(xintercept = c(-1, 1),
             linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") +
  theme_classic(base_size = 14) +
  labs(
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value",
    color = NULL,
    title = "OWA vs AM (F11)"
    ) +
  theme(plot.title = element_text(
    hjust = 0.5,
    size = 18,
    face = "bold"
      ))

volcano_df_F11

ggarrange(volcano_df_F0, volcano_df_F4, volcano_df_F11, nrow = 3)


##Make PCA with F0,4,11 gens with only OWA and AM####
vsd_AM_OWA <- vsd[, vsd$treatment %in% c("AM", "OWA")]

# Run PCA using all generations together
data <- plotPCA(
  vsd_AM_OWA,
  intgroup = c("treatment", "generation"),
  returnData = TRUE
)

# Calculate variance explained
percentVar <- round(100 * attr(data, "percentVar"))

# Plot all generations together

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
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  theme(legend.title = element_blank())
  

PCA_AM_OWA

##make upset plot####

res_OWAvsAM_F0 <- results(dds_F0, name="treatment_OWA_vs_AM", alpha=0.05)
res_OWAvsAM_F0 <- res_OWAvsAM_F0[order(res_OWAvsAM_F0$padj),]
res_OWAvsAM_F0 <- res_OWAvsAM_F0[!is.na(res_OWAvsAM_F0$padj),]
degs_OWAvsAM_F0 <- row.names(res_OWAvsAM_F0[res_OWAvsAM_F0$padj < 0.05,])

res_OWAvsAM_F4 <- results(dds_F4, name="treatment_OWA_vs_AM", alpha=0.05)
res_OWAvsAM_F4 <- res_OWAvsAM_F4[order(res_OWAvsAM_F4$padj),]
res_OWAvsAM_F4 <- res_OWAvsAM_F4[!is.na(res_OWAvsAM_F4$padj),]
degs_OWAvsAM_F4 <- row.names(res_OWAvsAM_F4[res_OWAvsAM_F4$padj < 0.05,])

res_OWAvsAM_F11 <- results(dds_F11, name="treatment_OWA_vs_AM", alpha=0.05)
res_OWAvsAM_F11 <- res_OWAvsAM_F11[order(res_OWAvsAM_F11$padj),]
res_OWAvsAM_F11 <- res_OWAvsAM_F11[!is.na(res_OWAvsAM_F11$padj),]
degs_OWAvsAM_F11 <- row.names(res_OWAvsAM_F11[res_OWAvsAM_F11$padj < 0.05,])

all_genes <- unique(c(
  degs_OWAvsAM_F0,
  degs_OWAvsAM_F4,
  degs_OWAvsAM_F11
))

upset_df <- data.frame(
  gene = all_genes,
  OWA = all_genes %in% degs_OWAvsAM_F0,
  OWA = all_genes %in% degs_OWAvsAM_F4,
  OWA = all_genes %in% degs_OWAvsAM_F11
)

head(upset_df)

deg.list <- list(
  F0 = degs_OWAvsAM_F0,
  F4 = degs_OWAvsAM_F4,
  F11 = degs_OWAvsAM_F11
)

upset(
  fromList(deg.list),
  order.by = "freq",
  mainbar.y.label = "Number of DEGs",
  sets.x.label = "Total DEGs"
)

##prettier
upset(
  fromList(deg.list),
  order.by = "freq",
  main.bar.color = "#CC3333",
  sets.bar.color = c("#CC3333", "#CC3333", "#CC3333"), # had to manually adjust the order
  mainbar.y.label = "Number of DEGs",
  sets.x.label = "Total DEGs"
)

##GO analysis!!####

##first make saved results files w trinity ids
#made for OWA vs AM already

#F0
res_OWAvsAM_F0.df <- as.data.frame(res_OWAvsAM_F0)
res_OWAvsAM_F0.df$fullID <- rownames(res_OWAvsAM_F0.df)

parts <- strsplit(res_OWAvsAM_F0.df$fullID, "::")

res_OWAvsAM_F0.df$shortID <- sapply(
  parts,
  function(x) paste(x[1:2], collapse="::")
)

write.csv(
  res_OWAvsAM_F0.df,
  "myresults/OWAvsAM_F0_results.csv",
  row.names = FALSE
)
#F4
res_OWAvsAM_F4.df <- as.data.frame(res_OWAvsAM_F4)
res_OWAvsAM_F4.df$fullID <- rownames(res_OWAvsAM_F4.df)

parts <- strsplit(res_OWAvsAM_F4.df$fullID, "::")

res_OWAvsAM_F4.df$shortID <- sapply(
  parts,
  function(x) paste(x[1:2], collapse="::")
)

write.csv(
  res_OWAvsAM_F4.df,
  "myresults/OWAvsAM_F4_results.csv",
  row.names = FALSE
)
#F11
res_OWAvsAM_F11.df <- as.data.frame(res_OWAvsAM_F11)
res_OWAvsAM_F11.df$fullID <- rownames(res_OWAvsAM_F11.df)

parts <- strsplit(res_OWAvsAM_F11.df$fullID, "::")

res_OWAvsAM_F11.df$shortID <- sapply(
  parts,
  function(x) paste(x[1:2], collapse="::")
)

write.csv(
  res_OWAvsAM_F11.df,
  "myresults/OWAvsAM_F11_results.csv",
  row.names = FALSE
)

#skip to here if already made csv's for GO####
##prep to run GO analyses, make sure I'm in mydata in working directory
setwd("/gpfs1/home/m/c/mcampuza/Projects/eco_genomics_2026/transcriptomics/mydata")
mappingFile <- "trinotate_annotation_GOblastx_forTopGO.txt"

geneID2GO <- readMappings(
  file = mappingFile
)

cat("Genes in GO mapping:",
    length(geneID2GO),
    "\n")

##Filter GO terms by size, show most relevant ones

allGO <- table(unlist(geneID2GO))

keepTerms <- names(allGO)[
  allGO >= 5 &
    allGO <= 500
]

geneID2GO.filtered <- lapply(
  geneID2GO,
  function(x) intersect(x, keepTerms)
)

geneID2GO.filtered <- geneID2GO.filtered[
  lengths(geneID2GO.filtered) > 0
]

##Create function to run TopGO contrast

run_topGO_contrast <- function(
    infile,
    outfile,
    ontology = "BP",
    padj.cutoff = 0.05){
  
  deseq <- read.csv(
    infile,
    stringsAsFactors = FALSE
  )
  
  deseq <- subset(
    deseq,
    !is.na(padj)
  )
  
  geneList <- factor(
    as.integer(deseq$padj < padj.cutoff)
  )
  
  names(geneList) <- deseq$shortID
  
  cat("\nGenes tested:",
      length(geneList))
  
  cat("\nSignificant genes:",
      sum(geneList == 1),
      "\n")
  
  GOdata <- new(
    "topGOdata",
    ontology = ontology,
    allGenes = geneList,
    geneSelectionFun = function(x) x == 1,
    annot = annFUN.gene2GO,
    gene2GO = geneID2GO.filtered
  )
  
  resultWeight <- runTest(
    GOdata,
    algorithm = "weight01",
    statistic = "fisher"
  )
  
  GOresults <- GenTable(
    GOdata,
    weightFisher = resultWeight,
    orderBy = "weightFisher",
    topNodes = 100
  )
  
  write.csv(
    GOresults,
    outfile,
    row.names = FALSE
  )
  
  return(GOresults)
}

##Run TopGO for different generations remember to set wd back to transcriptomics

setwd("/gpfs1/home/m/c/mcampuza/Projects/eco_genomics_2026/transcriptomics")

GO_OWA_F0 <- run_topGO_contrast(
  "myresults/OWAvsAM_F0_results.csv",
  "myresults/GO_OWAvsAM_F0_BP.csv"
)

GO_OWA_F4 <- run_topGO_contrast(
  "myresults/OWAvsAM_F4_results.csv",
  "myresults/GO_OWAvsAM_F4_BP.csv"
)

GO_OWA_F11 <- run_topGO_contrast(
  "myresults/OWAvsAM_F11_results.csv",
  "myresults/GO_OWAvsAM_F11_BP.csv"
)

# Read TopGO results
##F0
go_F0 <- read.csv(
  "myresults/GO_OWAvsAM_F0_BP.csv",
  stringsAsFactors = FALSE
)

# Convert p-values to numeric
go_F0$weightFisher <- gsub("^<\\s*", "", go_F0$weightFisher)

go_F0$weightFisher <- as.numeric(go_F0$weightFisher)

# Convert counts to numeric
go_F0$Significant <- as.numeric(go_F0$Significant)

# Create -log10(p)
go_F0$minusLogP <- -log10(go_F0$weightFisher)

# Keep top 10 GO terms
go_F0_top10 <- go_F0 %>%
  arrange(weightFisher) %>%
  slice(1:10)

# Order terms for plotting
go_F0_top10$Term <- factor(
  go_F0_top10$Term,
  levels = rev(go_F0_top10$Term)
)

# Bubble plot
GO_bubbleplot_F0_OWAvsAM <- ggplot(
  go_F0_top10,
  aes(
    x = minusLogP,
    y = Term
  )
) +
  geom_point(
    aes(
      size = Significant/Annotated,
      color = minusLogP
    )
  ) +
  scale_color_viridis_c() +
  theme_bw(base_size = 12) +
  labs(
    title = "Top GO Terms: OWA vs AM (F0)",
    x = expression(-logp),
    y = "GO Term",
    color = expression(-logp),
    size = "Proportion Significant\nGenes"
  )
GO_bubbleplot_F0_OWAvsAM

#F4
go_F4 <- read.csv(
  "myresults/GO_OWAvsAM_F4_BP.csv",
  stringsAsFactors = FALSE
)

# Convert p-values to numeric
go_F4$weightFisher <- gsub("^<\\s*", "", go_F4$weightFisher)

go_F4$weightFisher <- as.numeric(go_F4$weightFisher)

# Convert counts to numeric
go_F4$Significant <- as.numeric(go_F4$Significant)

# Create -log10(p)
go_F4$minusLogP <- -log10(go_F4$weightFisher)

# Keep top 10 GO terms
go_F4_top10 <- go_F4 %>%
  arrange(weightFisher) %>%
  slice(1:10)

# Order terms for plotting
go_F4_top10$Term <- factor(
  go_F4_top10$Term,
  levels = rev(go_F4_top10$Term)
)

# Bubble plot
GO_bubbleplot_F4_OWAvsAM <- ggplot(
  go_F4_top10,
  aes(
    x = minusLogP,
    y = Term
  )
) +
  geom_point(
    aes(
      size = Significant/Annotated,
      color = minusLogP
    )
  ) +
  scale_color_viridis_c() +
  theme_bw(base_size = 12) +
  labs(
    title = "Top GO Terms: OWA vs AM (F4)",
    x = expression(-logp),
    y = "GO Term",
    color = expression(-logp),
    size = "Proportion Significant\nGenes"
  )
GO_bubbleplot_F4_OWAvsAM

#F11
go_F11 <- read.csv(
  "myresults/GO_OWAvsAM_F11_BP.csv",
  stringsAsFactors = FALSE
)

# Convert p-values to numeric
go_F11$weightFisher <- gsub("^<\\s*", "", go_F11$weightFisher)

go_F11$weightFisher <- as.numeric(go_F11$weightFisher)

# Convert counts to numeric
go_F11$Significant <- as.numeric(go_F11$Significant)

# Create -log10(p)
go_F11$minusLogP <- -log10(go_F11$weightFisher)

# Keep top 10 GO terms
go_F11_top10 <- go_F11 %>%
  arrange(weightFisher) %>%
  slice(1:10)

# Order terms for plotting
go_F11_top10$Term <- factor(
  go_F11_top10$Term,
  levels = rev(go_F11_top10$Term)
)

# Bubble plot
GO_bubbleplot_F11_OWAvsAM <- ggplot(
  go_F11_top10,
  aes(
    x = minusLogP,
    y = Term
  )
) +
  geom_point(
    aes(
      size = Significant/Annotated,
      color = minusLogP
    )
  ) +
  scale_color_viridis_c() +
  theme_bw(base_size = 12) +
  labs(
    title = "Top GO Terms: OWA vs AM (F11)",
    x = expression(-logp),
    y = "GO Term",
    color = expression(-logp),
    size = "Proportion Significant\nGenes"
  )
GO_bubbleplot_F11_OWAvsAM

#put together
ggarrange(GO_bubbleplot_F0_OWAvsAM,GO_bubbleplot_F4_OWAvsAM,GO_bubbleplot_F11_OWAvsAM)

