library(topGO)
library(GO.db)
library(Rgraphviz)

# Please prepare the gene GO map file before running this script
geneID2GO <- readMappings('annotation/new_annotation_gomap.txt')

geneNames <- names(geneID2GO)

gene_interest <- geneNames[geneNames %in% filter(comp_st_26vs20_sig_no_st20b, logFC>0)$genes]

# define whether the gene present in the interest gene list
geneList <- factor(as.integer(geneNames %in% gene_interest))
names(geneList) <- geneNames

GOdata <- new("topGOdata", ontology = "MF", allGenes = geneList,
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

write_tsv(all_res, 'enrichment_analysis/st_26vs20_no_st20b.up.topgo.MF.txt')

all_res %>%
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
  ggtitle("st_26vs20_no_st20b up MF") + 
  scale_fill_distiller(palette = "Spectral")
