library(edgeR)
library(ggfortify)

# Load the expression table
pm_raw <- read.table(file="./CPm.count", sep="\t", row.names=1, header = T)

##################################################################
#   DE analysis using edgeR
##################################################################

# Make the DGElist object for the data
pm <- DGEList(counts = pm_raw
              , group = c(rep("12day", 2), rep("18day", 2), "24day")
              , genes = rownames(pm_raw)
)

# Estimates normalisation factors and dispersion
pm <- calcNormFactors(pm)  
y <- estimateCommonDisp(pm)
y <- estimateTagwiseDisp(y)

# Test the difference expression between groups
et <- exactTest(y,pair=c("12day","18day"))
et_12_24 <- exactTest(y,pair=c("12day","24day"))
et_18_24 <- exactTest(y,pair=c("18day","24day"))

# Write the DE analysis results
write.table(topTags(et, n = "all"), "./pm_12day_vs_18day_de.txt",sep="\t", quote = F, row.names = F)
write.table(topTags(et_12_24, n = "all"), "./pm_12day_vs_24day_de.txt",sep="\t", quote = F, row.names = F)
write.table(topTags(et_18_24, n = "all"), "./pm_18day_vs_24day_de.txt",sep="\t", quote = F, row.names = F)

###################################################################
#    PCA visualization
###################################################################
# Calculate the PCs based on the pseudo counts i.e. normalized count
pm.pseudocount.pca <- prcomp(t(y$pseudo.counts))
pm.rawcount.pca <- prcomp(t(pm_raw))
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
ggsave(plot = p, "pm.pseudocount.pca.pdf", width = 8, height = 7)


