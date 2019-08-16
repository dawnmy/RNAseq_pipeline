library(edgeR)
library(ggfortify)
library(tidyverse)
library(data.table)

# Load the expression table
pm_ko_raw <- read.table(file="./CPm.diamond.ko.count", sep="\t", row.names=1, header = T)
pm_ko_wo_noko <- pm_ko_raw[rownames(pm_ko_raw)!="noko",]
row.names(pm_ko_wo_noko) <- gsub("^ko:", "", row.names(pm_ko_wo_noko))
ko_func <- fread("../data/annotation/ko_annotation.txt", sep="@", header = F)
colnames(ko_func) <- c("genes", "function")

ko_func_uniq <- as.data.frame(ko_func[!duplicated(ko_func$genes),])
row.names(ko_func_uniq) <- ko_func_uniq$genes

pm_ko_func <- ko_func_uniq[row.names(pm_ko_wo_noko),]

##################################################################
#   DE analysis using edgeR
##################################################################

# Make the DGElist object for the data
pm <- DGEList(counts = pm_ko_wo_noko
              , group = c(rep("12day", 2), rep("18day", 2), "24day")
              , genes = rownames(pm_ko_wo_noko)
)

# Estimates normalisation factors and dispersion
pm <- calcNormFactors(pm)  
y <- estimateCommonDisp(pm)
y <- estimateTagwiseDisp(y)

# Test the difference expression between groups
et_18_12 <- exactTest(y,pair=c("18day","12day"))
et_24_12 <- exactTest(y,pair=c("24day","12day"))
et_24_18 <- exactTest(y,pair=c("24day","18day"))

# Add KEGG function for echo KO gene
add_func <- function(de){
  de_tab <- as.data.frame(topTags(de, n = "all"))
  de_func_tab <- cbind(de_tab, pm_ko_func[row.names(de_tab), ])
  return(de_func_tab)
}

# Write the DE analysis results
write.table(add_func(et_18_12), "./pm_18day_vs_12day_ko_de.txt",sep="\t", quote = F, row.names = F)
write.table(add_func(et_24_12), "./pm_24day_vs_12day_ko_de.txt",sep="\t", quote = F, row.names = F)
write.table(add_func(et_24_18), "./pm_24day_vs_18day_ko_de.txt",sep="\t", quote = F, row.names = F)

###################################################################
#    PCA visualization
###################################################################
# Calculate the PCs based on the pseudo counts i.e. normalized count
pm.pseudocount.pca <- prcomp(t(y$pseudo.counts))
pca_sample <- rownames(pm$samples)

# Visualize the PCA
p <- autoplot(pm.pseudocount.pca, 
              data=pm$samples, 
              colour='group', 
              size=5) + 
  theme_bw(base_size = 15) + 
  guides(shape = guide_legend(override.aes = list(size = 5)), 
         colour = guide_legend(override.aes = list(size = 5)),
         size = F)

# Save the figure
ggsave(plot = p, "pm.pseudocount.ko.pca.pdf", width = 8, height = 7)



###################################################################
#    KEGG enrichment analysis
###################################################################
library(clusterProfiler)
library(cowplot)
# The mothod can be "kover", "kgse", "kmover", "kmgse"
keggenrich_genes <- function(de, fc="up", method="kover"){
  if (fc == "up"){
    genes <- as.data.frame(topTags(de, n = "all")) %>% 
      filter(FDR<=0.05, logFC >= 1) %>% 
      .$genes %>% 
      gsub("^ko:", "", .)
  }
  else{
    genes = genes <- as.data.frame(topTags(de, n = "all")) %>% 
      filter(FDR<=0.05, logFC <= -1) %>% 
      .$genes %>% 
      gsub("^ko:", "", .)
  }
  # KEGG over representation analysis, Fisher's extact test
  if (method == "kover"){
    enrich <- enrichKEGG(gene         = genes, 
                         organism     = 'ko', 
                         pvalueCutoff = 0.05)
  }
  # KEGG enrichment analysis, KS test
  else if (method == "kgse") {
    enrich <- gseKEGG(geneList     = genes,
                      organism     = 'ko',
                      nPerm        = 1000,
                      minGSSize    = 25,
                      pvalueCutoff = 0.05,
                      verbose      = FALSE)
  } 
  # KEGG module over representation analysis
  else if (method == "kmover") {
    
    enrich <- enrichMKEGG(gene         = genes,
                          organism     = 'ko',
                          pvalueCutoff = 0.05)
  }
  # KEGG module enrichment analysis
  else{
    enrich <- gseMKEGG(geneList     = genes,
                       species      = 'ko',
                       nPerm        = 1000,
                       minGSSize    = 25,
                       pvalueCutoff = 0.05,
                       verbose      = FALSE)
  }
}

## Perform KEGG enrichment analysis using Fisher's exact test
kegg_enrich_up_24_12 <- keggenrich_genes(et_24_12, fc="up", method="kover")
kegg_enrich_down_24_12 <- keggenrich_genes(et_24_12, fc="down", method="kover")
p1 <- dotplot(kegg_enrich_up_24_12, showCategory=30) + ggtitle("KEGG enrich for up-regulated in 24day compared with 12day")
p2 <- dotplot(kegg_enrich_down_24_12, showCategory=30) + ggtitle("KEGG enrich for down-regulated in 24day compared with 12day")

ggsave("./kegg_enrich_24_12.pdf", 
       plot=plot_grid(p1, p2, ncol=1), 
       width = 8, 
       height = 10)

kegg_enrich_up_18_12 <- keggenrich_genes(et_18_12, fc="up", method="kover")
kegg_enrich_down_18_12 <- keggenrich_genes(et_18_12, fc="down", method="kover")
p3 <- dotplot(kegg_enrich_up_18_12, showCategory=30) + 
  ggtitle("KEGG enrich for up-regulated in 18day compared with 12day")
p4 <- dotplot(kegg_enrich_down_18_12, showCategory=30) + 
  ggtitle("KEGG enrich for down regulated in 18day compared with 12day")

ggsave("./kegg_enrich_18_12.pdf", 
       plot=plot_grid(p3, p4, ncol=1), 
       width = 8, 
       height = 10)

kegg_enrich_up_24_18 <- keggenrich_genes(et_24_18, fc="up", method="kover")
kegg_enrich_down_24_18 <- keggenrich_genes(et_24_18, fc="down", method="kover")
p5 <- dotplot(kegg_enrich_up_24_18, showCategory=30) + 
  ggtitle("KEGG enrich for up regulated in 24day compared with 18day")
p6 <- dotplot(kegg_enrich_down_24_18, showCategory=30) + 
  ggtitle("KEGG enrich for down regulated in 24day compared with 18day")

ggsave("./kegg_enrich_24_18.pdf", 
       plot=plot_grid(p5, p6, ncol=1), 
       width = 8, 
       height = 10)
