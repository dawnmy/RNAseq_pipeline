library(ggfortify)
library(tidyverse)

# Run this after DE script
pm.normcount.pca <- prcomp(t(cpm(y))
pm.rel.pca <- prcomp(t(pm_rel))
pca_sample <- rownames(y$samples)
metadata <- y$samples
metadata$sample <- rownames(y$samples)
metadata$rep <- substr(metadata$sample, 7,8)

# Visualize the PCA
pca_plot <- function(pc.tab, metadata){
  p <- autoplot(pc.tab, 
                data=metadata, 
                colour='group', 
                shape='rep',
                size=5) + 
    theme_bw(base_size = 15) + 
    guides(shape = guide_legend(override.aes = list(size = 5)), 
           colour = guide_legend(override.aes = list(size = 5)))  
}

pm.normalized.pcaplot <- pca_plot(pm.normcount.pca, metadata) + 
  ggtitle("PCA based on gene exp (normalized)")
pm.rel.pcaplot <- pca_plot(pm.rel.pca, metadata) + 
  ggtitle("PCA based on gene exp (relative abundance)")


# Save the figure
ggsave(plot=plot_grid(pm.normalized.pcaplot, pm.rel.pcaplot, ncol=2), 
       "pm.temperature.pseudocount.rel.pca.pdf", 
       width = 12, 
       height = 5)
