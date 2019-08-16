library(edgeR)
library(ggfortify)
library(tidyverse)
library(cowplot)

# Load the expression table
pm_raw <- read.table(file="./CPm.count", sep="\t", row.names=1, header = T)
pm_ <- pm_raw[rownames(pm_raw)!="*",]
pm_rel <- t(100*t(pm_)/colSums(pm_))


##################################################################
#   DE analysis using edgeR
##################################################################

# Make the DGElist object for the data
# Define the groups
pm <- DGEList(counts = pm_
              , group = c(rep("12day", 2), rep("18day", 2), "24day")
              , genes = rownames(pm_)
)

# Estimates normalisation factors and dispersion
pm <- calcNormFactors(pm)  
y <- estimateCommonDisp(pm)
y <- estimateTagwiseDisp(y)

# Test the difference expression between groups
et_18_12 <- exactTest(y,pair=c("18day","12day"))
et_24_12 <- exactTest(y,pair=c("24day","12day"))
et_24_18 <- exactTest(y,pair=c("24day","18day"))


# Write the DE analysis results
write.table(add_func(et_18_12), "./pm_18day_vs_12day_de.txt",sep="\t", quote = F, row.names = F)
write.table(add_func(et_24_12), "./pm_24day_vs_12day_de.txt",sep="\t", quote = F, row.names = F)
write.table(add_func(et_24_18), "./pm_24day_vs_18day_de.txt",sep="\t", quote = F, row.names = F)

###################################################################
#    PCA visualization
###################################################################
# Calculate the PCs based on the pseudo counts i.e. normalized count
pm.pseudocount.pca <- prcomp(t(y$pseudo.counts))
pm.rawcount.pca <- prcomp(t(pm_raw))
pm.rel.pca <- prcomp(t(pm_rel))
pca_sample <- rownames(pm$samples)

# Visualize the PCA
pca_plot <- function(pc.tab, metadata){
  p <- autoplot(pc.tab, 
           data=metadata, 
           colour='group', 
           size=5) + 
    theme_bw(base_size = 15) + 
    guides(shape = guide_legend(override.aes = list(size = 5)), 
           colour = guide_legend(override.aes = list(size = 5)),
           size = F)  
}

pm.normalized.pcaplot <- pca_plot(pm.pseudocount.pca, pm$samples) + 
  ggtitle("PCA based on gene exp (normalized)")
pm.rel.pcaplot <- pca_plot(pm.rel.pca, pm$samples) + 
  ggtitle("PCA based on gene exp (relative abundance)")


# Save the figure
ggsave(plot=plot_grid(pm.normalized.pcaplot, pm.rel.pcaplot, ncol=1), "pm.pseudocount.rel.pca.pdf", width = 7, height = 12)

