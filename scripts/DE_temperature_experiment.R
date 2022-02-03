library(edgeR)
library(ggfortify)
library(tidyverse)
library(cowplot)


# Load the expression table
pm_raw <- read.table(file="./count/pminimum.gene.feacturecount.allsamples.count.txt", sep="\t", 
                     row.names=1, 
                     header = T,
                     check.names = F) 

# Calculate the realative abundance
pm_rel <- t(100*t(pm_raw)/colSums(pm_raw))


##################################################################
#   DE analysis using edgeR
##################################################################

# Make the DGElist object for the data
# Define the groups
y <- DGEList(counts = pm_raw, 
             group = c(rep("20_ex", 3), rep("20_st", 3), rep("26_ex", 3),
                          rep("26_st", 3), rep("30_ex", 3), rep("30_st", 3)), 
             genes = rownames(pm_raw)
)


# define design matrix
temperature <- c(rep('20', 6),rep('26', 6),rep('30', 6))
phase <- c(rep('ex', 3), rep('st', 3), 
           rep('ex', 3), rep('st', 3),
           rep('ex', 3),rep('st', 3))

group <- factor(paste(phase, temperature, sep='_'))

design <- model.matrix(~0+group)

colnames(design) <- levels(group)

# Estimates normalisation factors and dispersion
y <- calcNormFactors(y)
y <- estimateDisp(y, design, robust=TRUE)


# y$common.dispersion
plotMDS(y, col=rep(1:6, each=3))
plotBCV(y)

# fit the design matrix
fit <- glmQLFit(y, design)
fit <- glmQLFit(y, design, robust=TRUE)
plotQLDisp(fit)

con_all <- makeContrasts(
  comp_ex_26vs20 = ex_26 - ex_20,
  comp_ex_30vs20 = ex_30 - ex_20,
  comp_ex_30vs26 = ex_30 - ex_26,
  comp_st_26vs20 = st_26 - st_20,
  comp_st_30vs20 = st_30 - st_20,
  comp_st_30vs26 = st_30 - st_26,
  comp_st_ex_20 = st_20 - ex_20, 
  comp_st_ex_26 = st_26 - ex_26,
  comp_st_ex_30 = st_30 - ex_30,
  levels=design
)

# pairwise test
# comp_ex_26vs20
qlf_comp_ex_26vs20 <- glmQLFTest(fit, contrast=con_all[,'comp_ex_26vs20'])
comp_ex_26vs20_sig <- topTags(qlf_comp_ex_26vs20, n='all') %>% 
  as.data.frame () %>%
  filter(FDR<=0.001, logCPM>=0, abs(logFC)>=2)

summary(decideTests(qlf_comp_ex_26vs20))
plotMD(qlf_comp_ex_26vs20)
abline(h=c(-2,2), col="purple")
abline(v=c(0), col="green")

# comp_ex_30vs20
qlf_comp_ex_30vs20 <- glmQLFTest(fit, contrast=con_all[,'comp_ex_30vs20'])
comp_ex_30vs20_sig <- topTags(qlf_comp_ex_30vs20, n='all') %>% 
  as.data.frame () %>%
  filter(FDR<=0.001, logCPM>=0, abs(logFC)>=2)

summary(decideTests(qlf_comp_ex_30vs20))
plotMD(qlf_comp_ex_30vs20)
abline(h=c(-2,2), col="purple")
abline(v=c(0), col="green")

# comp_st_26vs20
qlf_comp_st_26vs20 <- glmQLFTest(fit, contrast=con_all[,'comp_st_26vs20'])
comp_st_26vs20_sig <- topTags(qlf_comp_st_26vs20, n='all') %>% 
  as.data.frame () %>%
  filter(FDR<=0.001, logCPM>=0, abs(logFC)>=2)

summary(decideTests(qlf_comp_st_26vs20))
plotMD(qlf_comp_st_26vs20)
abline(h=c(-2,2), col="purple")
abline(v=c(0), col="green")

# comp_st_30vs20
qlf_comp_st_30vs20 <- glmQLFTest(fit, contrast=con_all[,'comp_st_30vs20'])
comp_st_30vs20_sig <- topTags(qlf_comp_st_30vs20, n='all') %>% 
  as.data.frame () %>%
  filter(FDR<=0.001, logCPM>=0, abs(logFC)>=2)

summary(decideTests(qlf_comp_st_30vs20))
plotMD(qlf_comp_st_30vs20)
abline(h=c(-2,2), col="purple")
abline(v=c(0), col="green")


# comp_ex_30vs26
qlf_comp_ex_30vs26 <- glmQLFTest(fit, contrast=con_all[,'comp_ex_30vs26'])
comp_ex_30vs26_sig <- topTags(qlf_comp_ex_30vs26, n='all') %>% 
  as.data.frame () %>%
  filter(FDR<=0.001, logCPM>=0, abs(logFC)>=2)

summary(decideTests(qlf_comp_ex_30vs26))
plotMD(qlf_comp_ex_30vs26)
abline(h=c(-2,2), col="purple")
abline(v=c(0), col="green")


# comp_st_30vs26
qlf_comp_st_30vs26 <- glmQLFTest(fit, contrast=con_all[,'comp_st_30vs26'])
comp_st_30vs26_sig <- topTags(qlf_comp_st_30vs26, n='all') %>% 
  as.data.frame () %>%
  filter(FDR<=0.001, logCPM>=0, abs(logFC)>=2)

summary(decideTests(qlf_comp_st_30vs26))
plotMD(qlf_comp_st_30vs26)
abline(h=c(-2,2), col="purple")
abline(v=c(0), col="green")

# compare same temperature different phase
# comp_st_ex_20
qlf_comp_st_ex_20 <- glmQLFTest(fit, contrast=con_all[,'comp_st_ex_20'])
comp_st_ex_20_sig <- topTags(qlf_comp_st_ex_20, n='all') %>% 
  as.data.frame () %>%
  filter(FDR<=0.001, logCPM>=0, abs(logFC)>=2)

summary(decideTests(qlf_comp_st_ex_20))
plotMD(qlf_comp_st_ex_20)
abline(h=c(-2,2), col="purple")
abline(v=c(0), col="green")


# comp_st_ex_26
qlf_comp_st_ex_26 <- glmQLFTest(fit, contrast=con_all[,'comp_st_ex_26'])
comp_st_ex_26_sig <- topTags(qlf_comp_st_ex_26, n='all') %>% 
  as.data.frame () %>%
  filter(FDR<=0.001, logCPM>=0, abs(logFC)>=2)

summary(decideTests(qlf_comp_st_ex_26))
plotMD(qlf_comp_st_ex_26)
abline(h=c(-2,2), col="purple")
abline(v=c(0), col="green")


# comp_st_ex_30
qlf_comp_st_ex_30 <- glmQLFTest(fit, contrast=con_all[,'comp_st_ex_30'])
comp_st_ex_30_sig <- topTags(qlf_comp_st_ex_30, n='all') %>% 
  as.data.frame () %>%
  filter(FDR<=0.001, logCPM>=0, abs(logFC)>=2)

summary(decideTests(qlf_comp_st_ex_30))
plotMD(qlf_comp_st_ex_30)
abline(h=c(-2,2), col="purple")
abline(v=c(0), col="green")


# VennDiagram

library(VennDiagram)

# Prepare a palette of 3 colors with R colorbrewer:
library(RColorBrewer)
myCol <- brewer.pal(4, "Pastel2")

venn.diagram(
  x = list(comp_st_26vs20_sig$genes, 
           comp_st_30vs26_sig$genes, 
           comp_st_30vs20_sig$genes),
  category.names = c("st_26vs20" , "st_30vs26" , "st_30vs20"),
  filename = 'analyses/comp_st_sig_venn_diagramm.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 1280 , 
  width = 1280 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol[1:3],
  
  # Numbers
  cex = .9,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.9,
  cat.fontfamily = "sans"
)


venn.diagram(
  x = list(comp_st_26vs20_sig$genes, 
           comp_st_30vs26_sig$genes, 
           comp_ex_26vs20_sig$genes, 
           comp_ex_30vs26_sig$genes),
  category.names = c("st_26vs20" , "st_30vs26" , "ex_26vs20", "ex_30vs26"),
  filename = 'analyses/comp_vs26_sig_venn_diagramm.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 1280 , 
  width = 1480 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = .9,
  fontface = "bold",
  fontfamily = "sans",

  # Set names
  cat.cex = 0.9,
  cat.fontfamily = "sans",
)

# same temperature different phase
venn.diagram(
  x = list(comp_st_ex_20_sig$genes, 
           comp_st_ex_26_sig$genes, 
           comp_st_ex_30_sig$genes),
  category.names = c("st_ex_20" , "st_ex_26" , "st_ex_30"),
  filename = 'analyses/comp_st_vs_ex_sig_venn_diagramm.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 1280 , 
  width = 1280 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol[1:3],
  
  # Numbers
  cex = .9,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.9,
  cat.fontfamily = "sans",
)

