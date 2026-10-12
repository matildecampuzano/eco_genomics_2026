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


##prep to run GO analyses, make sure I'm in mydata in working directory

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
