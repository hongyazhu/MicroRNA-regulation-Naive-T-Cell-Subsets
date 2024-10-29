# PCA-related analysis for decay rates of naive T cell datasets
# Supp Fig 2G-I

decay_rates = read.table('inputs/exonIntron_matrix/decay_rates.txt')
sample_info = read.table('inputs/exonIntron_matrix/sample_info.txt',
                         sep = '\t')

decay_rates_10exon = read.table('inputs/exonIntron_matrix/decay_rates_match_noNA_noInf_noall0_10exon.txt')

decay_rates_pseudocount0_noNA <- na.omit(decay_rates)
decay_rates_pseudocount0_noNA_noInf <- decay_rates_pseudocount0_noNA[!is.infinite(rowSums(decay_rates_pseudocount0_noNA)),] # 6529 genes, mirna_gene_naiveOnly_moresamples220415/inputs/exonIntron_matrix/decay_rates_match_noNA_noInf.txt also 6529 genes
decay_rates_pseudocount0_noNA_noInf_10exon = decay_rates_pseudocount0_noNA_noInf[intersect(rownames(decay_rates_pseudocount0_noNA_noInf), rownames(decay_rates_10exon)),]

colSums(decay_rates_pseudocount0_noNA_noInf)

table(sample_info$Sample == colnames(decay_rates_pseudocount0_noNA_noInf))
cnames_check = gsub('-', '.', sample_info$Sample)
cnames_check[29:36] = paste0('X', cnames_check[29:36])
table(colnames(decay_rates_pseudocount0_noNA_noInf) == cnames_check)

rownames(sample_info) = colnames(decay_rates_pseudocount0_noNA_noInf)
# sample_info[sample_info == "Wildtype"] <- "WT"
# sample_info[sample_info == 'Cre+'] = 'CrePos'
# sample_info[sample_info == 'Cre-'] = 'CreNeg'

library(ggfortify)

decay_rates_pseudocount0_noNA_noInf_10exon_scaled_t = data.frame(scale(t(decay_rates_pseudocount0_noNA_noInf_10exon)))
pca<-prcomp(decay_rates_pseudocount0_noNA_noInf_10exon_scaled_t)
rownames(pca$x) =colnames(decay_rates_pseudocount0_noNA_noInf_10exon)
autoplot(pca, label = T, label.size = 3.5) 
# PC1: 26.8%
# PC2: 10.91%
autoplot(pca, label = T, label.size = 3.5, y = 3) 
# PC3: 8.41%
autoplot(pca, label = T, label.size = 3.5, y = 4) 
# PC4: 6.04%


pca_plot = as.data.frame(pca$x)
pca_plot$Sample = sample_info$Sample
pca_plot$Source = sample_info$Source
pca_plot$Age = sample_info$Age
pca_plot$CellType = sample_info$CellType
pca_plot$Organ = sample_info$Organ
pca_plot$Treatment = sample_info$Treatment
pca_plot$Clean = sample_info$Clean

library(RColorBrewer)
mycolors = c(brewer.pal(name="Set2", n = 8))
ggplot(pca_plot, aes(x = PC1, y= PC2)) +
  geom_point(size = 3, aes(colour = Source, shape = Age)) +
  theme_bw(base_size = 14) +
  xlab("PC1 (26.8%)") +
  ylab("PC2 (10.91%)") +
  labs(color='Source') +
  labs(shape='Age') +
  guides(color = guide_legend(order = 1),
         shape = guide_legend(order = 2)) +
  scale_color_manual(limits = c('Clean Project', 'Mir29 Project', 'Wang 2016', 'Wissink 2015', 
                                'Mir29Cre Project', 'CD8Dev Project'), # 'Clean Thymocytes', 
                     labels = c('Tabilas et al. (2022)', 'Yee Mon et al. (2021)', 'Wang et al. (2016)', 'Wissink et al. (2015)', 
                                'Mir29Cre Project', 'CD8Dev Project'), # 'Clean Thymocytes', 
                     values = mycolors) #+ scale_shape_discrete(labels = c('Bulk', 'VM', 'TN'))
ggsave("decayRates_PCA/PCA_vst_blindT_sourceAge.pdf", plot = last_plot(), device = "pdf",
       width = 6.5, height = 4, dpi = 300)

ggplot(pca_plot, aes(x = PC1, y= PC2)) +
  geom_point(size = 3, aes(colour = Age, shape = CellType)) +
  theme_bw(base_size = 14) +
  xlab("PC1 (26.8%)") +
  ylab("PC2 (10.91%)") +
  # labs(color='Source') +
  # labs(color='Cell type') +
  guides(color = guide_legend(order = 1),
         shape = guide_legend(order = 2)) +
  scale_color_manual(values = mycolors[3:4])
ggsave("decayRates_PCA/PCA_vst_blindT_ageCellType.pdf", plot = last_plot(), device = "pdf",
       width = 5.5, height = 4, dpi = 300)


library(PCAtools)
# repeat previous result
p0 <- pca(t(decay_rates_pseudocount0_noNA_noInf_10exon_scaled_t))
biplot(p0)
# add metadata
metadata = sample_info
metadata$Sample = NULL
metadata$Source = factor(metadata$Source, levels = unique(metadata$Source)) # randomly ordered
metadata$Age = factor(metadata$Age, levels = c("Adult", 'Neo')) 
#metadata$CellType <- replace(metadata$CellType, metadata$CellType == 'MP', 'VM') 
metadata$CellType = factor(metadata$CellType, levels = c("TN", 'Bulk','VM')) 
#metadata$Organ = factor(metadata$Organ, levels = c("Spleen", 'Thymus')) 
metadata$Treatment = factor(metadata$Treatment, levels = c("WT", 'CreNeg', 'iLin28b','CrePos')) 
metadata$Clean = factor(metadata$Clean, levels = c("Clean", 'Dirty')) 

metadata_tmp = metadata
colnames(metadata_tmp) = c('Source', 'Neonatal/Adult', 'VM/TN', 'Dirty/Clean') # 'Thymus/Spleen', 
metadata_tmp = metadata_tmp[,rev(c(3,2))]
p <- pca(t(decay_rates_pseudocount0_noNA_noInf_10exon_scaled_t), metadata = metadata_tmp)
biplot(p)
pdf("decayRates_PCA/PCA_eigencorplot.pdf", width = 7, height = 4)
eigencorplot(p,
             components = getComponents(p, 1:3),
             metavars = rev(c('VM/TN','Neonatal/Adult')), # 'Thymus/Spleen',
             col = c('white', 'cornsilk1', 'gold', 'forestgreen', 'darkgreen'),
             cexCorval = 1.2,
             fontCorval = 2,
             posLab = 'all',
             #rotLabX = 45,
             #scale = TRUE,
             main = bquote(PC ~ Pearson ~ r^2 ~ metadata ~ correlates),
             plotRsquared = TRUE,
             corFUN = 'pearson',
             corUSE = 'pairwise.complete.obs',
             corMultipleTestCorrection = 'BH',
             signifSymbols = c('****', '***', '**', '*', ''),
             signifCutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1))
dev.off()
