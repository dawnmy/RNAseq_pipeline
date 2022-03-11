# load gene length
gene_length <- read.delim('annotation/gene_length.txt', sep='\t', row.names = 1)

# make sample name list with ex, st order
ex_st_samples <- paste(c(rep("20_ex", 3), rep("26_ex", 3), rep("30_ex", 3), 
        rep("20_st", 3), rep("26_st", 3), rep("30_st", 3)), 
      rep(c('a','b','c'), 6), sep='_')

# ex and st samples
ex_samples <- ex_st_samples[1:9]
st_samples <- ex_st_samples[10:18]

# calculate centered logRPKM based on DEG object and gene length
rpkm <- rpkm(y, gene.length = gene_length[y$genes$genes,])
logrpkm <- log2(rpkm + 1)
centered_logrpkm <- logrpkm - rowMeans(logrpkm)

# compute the group RPKM
rpkm_by_group <- rpkmByGroup(y, gene.length = gene_length[y$genes$genes,])
# write.table(as.data.frame(rpkm_by_group), file='analyses/rpkm_by_group.txt', sep='\t', quote = F)

# centered logRPKM of all significantly DE genes in ex ANOVA comparison and st ANOVA comparison
ex_sig_centered_logrpkm <- centered_logrpkm[ex_anova_de_sig$genes, ex_samples]
st_sig_centered_logrpkm <- centered_logrpkm[st_anova_de_sig$genes, st_samples]

# hierachical clustering of the expression
ex_sig_clusters <- hclust(dist(ex_sig_centered_logrpkm))
st_sig_clusters <- hclust(dist(st_sig_centered_logrpkm))


ex_sig_cut4clusters <- cutree(ex_sig_clusters, 4)
st_sig_cut4clusters <- cutree(st_sig_clusters, 4)

ex_sig_cut8clusters <- cutree(ex_sig_clusters, 8)
ex_sig_cut15clusters <- cutree(ex_sig_clusters, 15)
st_sig_cut8clusters <- cutree(st_sig_clusters, 8)
# coolmap(logrpkm_logFC6_comp_st_26vs20_sig, margins=c(7,7), lhei=c(1,6), lwid=c(1,3))

# plot(ex_sig_clusters)

# make wide table for gene clusters
ex_sig_cut4clusters_wide <- left_join(rownames_to_column(as.data.frame(ex_sig_centered_logrpkm), var='genes'), 
          data.frame(clusters=ex_sig_cut4clusters, 
                     genes = names(ex_sig_cut4clusters))) %>%
  mutate(clusters=factor(paste0('SC', clusters),
                         levels=paste0('SC', sort(unique(clusters)))))

ex_sig_cut8clusters_wide <- left_join(rownames_to_column(as.data.frame(ex_sig_centered_logrpkm), var='genes'), 
                                      data.frame(clusters=ex_sig_cut8clusters, 
                                                 genes = names(ex_sig_cut8clusters))) %>%
  mutate(clusters=factor(paste0('SC', clusters),
                         levels=paste0('SC', sort(unique(clusters)))))

st_sig_cut4clusters_wide <- left_join(rownames_to_column(as.data.frame(st_sig_centered_logrpkm), var='genes'), 
                                      data.frame(clusters=st_sig_cut4clusters, 
                                                 genes = names(st_sig_cut4clusters))) %>%
  mutate(clusters=factor(paste0('SC', clusters),
                         levels=paste0('SC', sort(unique(clusters)))))

st_sig_cut8clusters_wide <- left_join(rownames_to_column(as.data.frame(st_sig_centered_logrpkm), var='genes'), 
                                      data.frame(clusters=st_sig_cut8clusters, 
                                                 genes = names(st_sig_cut8clusters))) %>%
  mutate(clusters=factor(paste0('SC', clusters),
         levels=paste0('SC', sort(unique(clusters)))))


# all in one function for clustering and output wide table of clusters for genes
make_de_sig_clustes_wide <- function(phase='ex', cuth = 6, method='complete'){
  
  if (phase=='ex'){
    de_sig_centered_logrpkm <- ex_sig_centered_logrpkm
  }
  else{
    de_sig_centered_logrpkm <- st_sig_centered_logrpkm
  }
 
 
  de_sig_clusters <- hclust(dist(de_sig_centered_logrpkm), method = method)
  de_sig_cutnclusters <- cutree(de_sig_clusters, h=cuth)
  out_clusters_wide <- left_join(rownames_to_column(as.data.frame(de_sig_centered_logrpkm), 
                                                           var='genes'), 
                                        data.frame(clusters=de_sig_cutnclusters, 
                                                   genes = names(de_sig_cutnclusters))) %>%
    mutate(clusters=factor(paste0('SC', clusters)))
  
  # Is did not work to define levels directly in the mutate above possibly due to a bug of dplyr
  # workaround here is to change the levels outside of the pipe on the dataframe object
  levels(out_clusters_wide$clusters) <- paste0('SC', 
                                               sort(as.numeric(gsub('SC', '', 
                                                                    unique(out_clusters_wide$clusters)))))
  return(out_clusters_wide)
  
}


# all in one function for visualizing the clusters with line plot
cluster_line_plot <- function(clusters_wide){
  mean_cluster_wide <- sapply(c("20", "26", "30"), 
         function(x) rowMeans(clusters_wide[, grep(paste0('^', x), names(clusters_wide))])) %>%
    as.data.frame()
  mean_cluster_wide$genes <- clusters_wide$genes
  mean_cluster_wide$clusters <- clusters_wide$clusters
  
  mean_clusters_long <- mean_cluster_wide %>%
    gather(key='samples', value='logRPKM', -c(genes, clusters)) %>%
    mutate(temperatures=substr(samples, 1,2))
  
  clusters_count <- group_by(mean_cluster_wide, clusters) %>%summarise(n=n()) %>%
    mutate(count=paste0('N = ', n))
  
  ggplot(mean_clusters_long, aes(temperatures, logRPKM, alpha=0.3, group=genes)) + 
    geom_point(alpha=0.3, size=0.8) +
    geom_line(color='#409BD2') +
    # geom_line(aes(group = genes), color='#409BD2') +
    stat_summary(aes(group=1), fun=mean, geom="line", colour="red", size=1.5) +
    facet_wrap(~clusters) +
    geom_text(data=clusters_count, aes(x=1.8, y=5, label=count), 
              colour="black", inherit.aes=FALSE, parse=FALSE) +
    theme_bw(base_size = 14) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          # Change axis line
          axis.line = element_line(colour = "black"),
          legend.position='none') +
    ylab('centered log2(RPKM+1)') +
    scale_x_discrete(expand = c(0.05, 0.05))
}

cluster_line_plot(ex_sig_cut4clusters_wide)
cluster_line_plot(ex_sig_cut8clusters_wide)
cluster_line_plot(st_sig_cut4clusters_wide)
cluster_line_plot(st_sig_cut8clusters_wide)

# directly plot did not work, have to generate plot object then save to file
plot1 <- cluster_line_plot(make_de_sig_clustes_wide(phase='ex', cuth=5, method = 'complete'))
plot2 <- cluster_line_plot(make_de_sig_clustes_wide(phase='st', cuth=4, method = 'complete'))
# cluster_temp <- make_de_sig_clustes_wide(phase='ex', cuth=2, method = 'complete')

ggsave('analyses/st_cluster_cuth4.pdf', plot=plot2)

## Heatmap of gene expression in RPKM
library(ComplexHeatmap)
library(RColorBrewer)

cluster_exp_heatmap <- function(clusters_wide, samples){
  clusters_rpkm <- rpkm[clusters_wide[order(clusters_wide$clusters),'genes'], samples]
  
  
  num_clusters <- length(unique(clusters_wide$clusters))
  column_ha = columnAnnotation(temperatures = substr(samples, 1, 2),
                               col = list(temperatures=setNames(brewer.pal(name = "Set1", n = 3), 
                                                                c('20', '26', '30')))
  )
  row_ha = rowAnnotation(
    clusters = clusters_wide[order(clusters_wide$clusters),'clusters'],
    col = list(clusters=setNames(brewer.pal(name = "Set3", n = num_clusters), 
                                 paste0('SC', 1:num_clusters)))
  )
  
  exp_heat <- Heatmap(t(scale(t(clusters_rpkm))), 
                      name = "z-scaled RPKM", 
                      top_annotation = column_ha, 
                      left_annotation = row_ha,
                      cluster_columns = F,
                      cluster_rows = F,
                      show_row_names = F,
                      col = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
                      # heatmap_legend_param = list(direction = "horizontal")
  )
  exp_heat
  
}

cluster_exp_heatmap(st_sig_cut8clusters_wide, st_samples)

pdf("analyses/allsample_aldex_vitamin_B12_gene_exp_genus_bar_new.pdf", height = 14.8, width = 14)
# draw(exp_heat, merge_legend = TRUE)
draw(exp_heat, merge_legend = TRUE, heatmap_legend_side = "bottom",
     annotation_legend_side = "bottom")
dev.off()



## Function enrich of clusters
enrich_on_cluster <- function(cluster_wide, cluster, ontology='MF'){
  gene_interest <- geneNames[geneNames %in% filter(cluster_wide, clusters==cluster)$genes]
  
  # define whether the gene present in the interest gene list
  geneList <- factor(as.integer(geneNames %in% gene_interest))
  names(geneList) <- geneNames
  
  GOdata <- new("topGOdata", ontology = ontology, allGenes = geneList,
                annot = annFUN.gene2GO, gene2GO = geneID2GO)
  
  # define test using the weight01 algorithm (default) with fisher
  weight_fisher_result <- runTest(GOdata, 
                                  algorithm='weight01', 
                                  statistic='fisher') 
  
  # generate a table of results: we can use the GenTable function to generate a summary table with the results from tests applied to the topGOdata object.
  allGO <- usedGO(GOdata)
  all_res <- GenTable(GOdata, 
                      weightFisher=weight_fisher_result, 
                      orderBy='weightFisher', 
                      topNodes=length(allGO),
                      numChar=1000)
  
  #performing BH correction on our p values
  pvalues <- unlist(lapply(all_res$weightFisher, 
                           function(x) ifelse(x=="< 1e-30", 1e-30, as.numeric(x))))
  
  all_res$p.adj <- round(p.adjust(pvalues,method="BH"),
                         digits = 8)
  
  all_res
}


cluster_enrich_go <- enrich_on_cluster(ex_sig_cut4clusters_wide, 'SC3', ontology='BP')

write_tsv(cluster_enrich_go, 'analyses/ex_sig_cut4clusters_wide.SC3.topgo.BP.txt')

cluster_enrich_go %>%
  filter(p.adj<=0.05, Significant>=5) %>%
  arrange(Significant/Annotated) %>%
  mutate(Term=factor(Term, levels=.$Term)) %>%
  ggplot(aes(Significant/Annotated, 
             Term, 
             size=Significant, 
             fill=p.adj)) +
  geom_point(pch=21) +
  theme_bw(base_size = 15) +
  ylab("") +
  ggtitle("ex_cut4clusters SC3 BP") + 
  scale_fill_distiller(palette = "Spectral")
