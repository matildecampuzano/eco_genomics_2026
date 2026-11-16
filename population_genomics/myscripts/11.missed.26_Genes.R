library(tidyverse)
library(topGO)
#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("topGO")

#outliers detected by XtX scan 
outlier_file <- "Lottia_bp_xtx_outliers.txt"  

##map gene to GO terms
map_file      <- "lottia_gene2GO.map"         
###region of genes 
snp2gene_file <- "lottia_snp2gene_genic.tsv" 

###read in GO terms 
gene2GO <- readMappings(map_file)
length(gene2GO)          # how many genes have at least one GO term?
gene2GO[1:2]             # what does an entry look like?

##creating the gene universe - all genes in experiment 
snp2gene <- read.delim(snp2gene_file, header = FALSE,
                       col.names = c("CHROM", "POS", "GeneId"))
universe <- unique(snp2gene$GeneId)
length(universe)

mean(universe %in% names(gene2GO))    # proportion of universe genes with GO terms

## Genes that contain at least one outlier SNP.
outliers   <- read.table(outlier_file, header = TRUE)
hits       <- merge(outliers[, c("CHROM", "POS")], snp2gene, by = c("CHROM", "POS"))
candidates <- unique(hits$GeneId)

##how many candiate genes do we have?
cat("Outlier SNPs:", nrow(outliers),
    "| outlier SNPs inside genes:", nrow(unique(hits[, c("CHROM", "POS")])),
    "| candidate genes:", length(candidates), "\n")

geneList <- factor(as.integer(universe %in% candidates), levels = c(0, 1))
names(geneList) <- universe
table(geneList)

## nodeSize = 10 ignores GO terms with fewer than 10 annotated genes (we can add a max here if we have very broad terms)
GOdata <- new("topGOdata",
              ontology = "BP", ##biological processes, what about other options?
              allGenes = geneList,
              annot    = annFUN.gene2GO,
              gene2GO  = gene2GO,
              nodeSize = 10)
GOdata   # read this summary: how many genes/candidates actually have GO terms?

## test for enrichment 
## Fisher's exact test asks: are candidates over-represented in this term compared with the universe?
## 'weight01' accounts for the GO hierarchy: if a specific child term explains the signal, its broad parent terms are down-weighted.
result <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
result
##results table 
## Annotated   = universe genes in the term
## Significant = candidate genes in the term
## Expected    = candidates expected by chance
res <- GenTable(GOdata, weight01 = result, orderBy = "weight01",
                topNodes = 50, numChar = 1000)
res$p <- as.numeric(sub("< ", "", res$weight01))   
head(res, 10)


##write our outliers
write.csv(res, "topGO_BP_xtx.csv", row.names = FALSE)

## Which candidate genes are in the top term?
top_term <- res$GO.ID[1]
intersect(genesInTerm(GOdata, top_term)[[1]], sigGenes(GOdata))

## GeneRatio = candidate genes in the term / all candidates with GO terms
plot_df <- subset(res, p < 0.05)
if (nrow(plot_df) == 0) {
  message("No terms with p < 0.05 - plotting the top 15 anyway (none significant)")
  plot_df <- head(res, 15)
}
plot_df <- head(plot_df, 20)
plot_df$GeneRatio <- plot_df$Significant / numSigGenes(GOdata)
plot_df$logP      <- -log10(plot_df$p)
## wrap long GO term names over two lines so the text doesn't dominate the plot
plot_df$Label <- vapply(plot_df$Term, function(t) paste(strwrap(t, width = 40), collapse = "\n"),
                        character(1))
plot_df$Label <- reorder(plot_df$Label, plot_df$GeneRatio)

p1 <- ggplot(plot_df, aes(x = Label, y = GeneRatio, size = Significant, fill = logP)) +
  geom_segment(aes(y = 0, yend = GeneRatio, x = Label, xend = Label),
               colour = "lightblue", linewidth = 1.5, alpha = 0.5) +
  geom_point(shape = 21, alpha = 0.8, colour = "black") +
  coord_flip() +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4)) +
  theme_bw(base_size = 12) +
  theme(axis.text.y = element_text(colour = "black", lineheight = 0.9),
        plot.title.position = "plot",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_fill_gradientn(colours = c("#2b83ba", "#fdae61", "#d7191c"),
                       name = "-log10(p value)") +
  scale_size(range = c(3, 10), name = "Count") +
  ylab("Gene Ratio") + xlab("")
p1

#By environmental factors!!!!####
library(tidyverse)
env <- outliers_BF

geneanno<-read.table("outlier_gene_effects.tsv")

head(geneanno)

names(geneanno) <- c("CHROM", "POS", "GeneID", "Annotation")

envanno<- env %>% left_join(geneanno, by = c("CHROM", "POS"))

envanno %>% filter(!is.na(Annotation)) ->envanno

ggplot(aes(x = Annotation, fill = Annotation), data = envanno) + 
  geom_bar(stat = "count")+
  theme_classic()+
  theme(legend.position = "none")+
  coord_flip()+
  scale_fill_brewer(palette = "Set3")+
  facet_wrap(~covariate)



## the same numbers as a table
envanno %>% count(covariate, Annotation) %>% arrange(covariate, desc(n))

envanno %>% filter(covariate == "ph_min") %>% arrange(desc(BF_dB)) %>% head()
envanno %>% filter(covariate == "thetao_max") %>% arrange(desc(BF_dB)) %>% head()
