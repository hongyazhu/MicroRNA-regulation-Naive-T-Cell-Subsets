
library(ggplot2)
outputdir = 'inputs/exonIntron_matrix/'


#############
### mir29 ###
#############

mir29_counts_exon = read.table("featureCounts_output/mir29EV_counts_exon.txt", header = T, row.names = 1)
exon_length = mir29_counts_exon$Length
mir29_counts_exon = mir29_counts_exon[,6:ncol(mir29_counts_exon)]
colnames(mir29_counts_exon) = c("NeoM29-3-VM", 
                                "Adult-1-TN", "Adult-1-VM", "NeoNC-1-TN", "NeoNC-1-VM", "NeoM29-1-TN", "NeoM29-1-VM", 
                                "Adult-2-TN", "Adult-2-VM", "NeoNC-2-TN", "NeoNC-2-VM", "NeoM29-2-TN", "NeoM29-2-VM",
                                "Adult-3-TN", "Adult-3-VM", "NeoNC-3-TN", "NeoNC-3-VM", "NeoM29-3-TN") 
mir29_counts_exon = mir29_counts_exon[, c("Adult-1-TN", "Adult-1-VM", "NeoNC-1-TN", "NeoNC-1-VM", 
                                          "Adult-2-TN", "Adult-2-VM", "NeoNC-2-TN", "NeoNC-2-VM", 
                                          "Adult-3-TN", "Adult-3-VM", "NeoNC-3-TN", "NeoNC-3-VM") ]
mir29_counts_exon$gene = rownames(mir29_counts_exon)


################
### erin2015 ###
################

erin2015_counts_exon = read.table("featureCounts_output/Wissink2015_counts_exon.txt", header = T, row.names = 1)
exon_length = erin2015_counts_exon$Length
erin2015_counts_exon = erin2015_counts_exon[,6:ncol(erin2015_counts_exon)]
colnames(erin2015_counts_exon) = gsub("X.home.hz543.data.Erin2015_RNAseq_all.Jan_2021.combine_fastq_to_sample.hisat2.", "", colnames(erin2015_counts_exon))
colnames(erin2015_counts_exon) = gsub(".bam", "", colnames(erin2015_counts_exon))
erin2015_counts_exon = erin2015_counts_exon[, c('adult_naive', 'neonate_naive')]
erin2015_counts_exon$gene = rownames(erin2015_counts_exon)


#################
### blood2016 ###
#################

blood2016_counts_exon = read.table("featureCounts_output/blood2016_counts_exon.txt", header = T, row.names = 1)
exon_length = blood2016_counts_exon$Length
blood2016_counts_exon = blood2016_counts_exon[,6:ncol(blood2016_counts_exon)]
colnames(blood2016_counts_exon) = sub(".bam", "", colnames(blood2016_counts_exon))
colnames(blood2016_counts_exon) = sub("X.home.hz543.data.2016Blood.neo.adu.April032021.hisat2.", "", colnames(blood2016_counts_exon))
blood2016_counts_exon = blood2016_counts_exon[, c('adult_spleen', 'neonatal_spleen')] #, 'adult_thymus', 'neonatal_thymus'
blood2016_counts_exon$gene = rownames(blood2016_counts_exon)


################
### cleanCd8 ###
################

cleanCd8_counts_exon = read.table("featureCounts_output/cleanCd8_counts_exon.txt", header = T, row.names = 1)
exon_length = cleanCd8_counts_exon$Length
cleanCd8_counts_exon = cleanCd8_counts_exon[,6:ncol(cleanCd8_counts_exon)]
coln_tmp = sub('.*CD', 'CD', colnames(cleanCd8_counts_exon))
coln_tmp = sub('_CKDL.*', '', coln_tmp)
targetFile = read.table('targetFile.txt', header = T)
targetFile$cd = sub('.ReadsPerGene.out.tab.rawCounts', '', targetFile$files)
colnames(cleanCd8_counts_exon) = targetFile[match(coln_tmp, targetFile$cd),]$label
cleanCd8_counts_exon = cleanCd8_counts_exon[, c('TN.C.4wts.1', 'TN.C.4wts.2', 'TN.C.4wts.3', 
                                                'VM.C.4wts.1', 'VM.C.4wts.2', 'VM.C.4wts.3',
                                                'TN.D.4wts.1', 'TN.D.4wts.2', 'TN.D.4wts.3',
                                                'VM.D.4wts.1', 'VM.D.4wts.2', 'VM.D.4wts.3')]
colnames(cleanCd8_counts_exon) = c("Clean_TN_rep1", "Clean_TN_rep2", "Clean_TN_rep3",  
                                   "Clean_VM_rep1", "Clean_VM_rep2", "Clean_VM_rep3",
                                   "Dirty_TN_rep1", "Dirty_TN_rep2", "Dirty_TN_rep3", 
                                   "Dirty_VM_rep1", "Dirty_VM_rep2", "Dirty_VM_rep3")
cleanCd8_counts_exon$gene = rownames(cleanCd8_counts_exon)


################
### mir29cre ###
################

mir29cre_counts_exon = read.table("featureCounts_output/mir29cre_counts_exon.txt", header = T, row.names = 1)
exon_length = mir29cre_counts_exon$Length
mir29cre_counts_exon = mir29cre_counts_exon[,6:ncol(mir29cre_counts_exon)]
colnames(mir29cre_counts_exon) = sub("X.home.hz543.data.1087.Ciaran.mir29KO.RNAseq.redo_220406.hisat2.", "", colnames(mir29cre_counts_exon))
colnames(mir29cre_counts_exon) = c('1.WT.TN', '1.WT.MP', '1.mCre.TN', '1.mCre.MP', '1.pCre.TN', '1.pCre.MP',
                                   '2.WT.TN', '2.WT.MP', '2.mCre.TN', '2.mCre.MP', '2.pCre.TN', '2.pCre.MP')
mir29cre_counts_exon = mir29cre_counts_exon[, c('1.WT.TN', '1.WT.MP', '1.mCre.TN', '1.mCre.MP',
                                                '2.WT.TN', '2.WT.MP', '2.mCre.TN', '2.mCre.MP')]
mir29cre_counts_exon$gene = rownames(mir29cre_counts_exon)


##############
### cd8dev ###
##############

library("readxl")
filenames <- read_excel('cd8dev_Filenames.xlsx', col_names = F)
filenames$DGnum = gsub('N.ReadsPerGene.out.tab.rawCounts', '', filenames$...3)
filenames$DGnum = gsub('c', '', filenames$DGnum)

cd8dev_counts_exon = read.table("featureCounts_output/cd8dev_counts_exon.txt", header = T, row.names = 1)
exon_length = cd8dev_counts_exon$Length
cd8dev_counts_exon = cd8dev_counts_exon[,6:ncol(cd8dev_counts_exon)]
colnames(cd8dev_counts_exon) = gsub('.*.DG', 'DG', colnames(cd8dev_counts_exon))
colnames(cd8dev_counts_exon) = gsub('N.*', '', colnames(cd8dev_counts_exon))
#table(as.data.frame(filenames[match(colnames(cd8dev_counts_exon),filenames$DGnum),'DGnum'])$DGnum == colnames(cd8dev_counts_exon))
colnames(cd8dev_counts_exon) = as.data.frame(filenames[match(colnames(cd8dev_counts_exon),filenames$DGnum),'...2'])$'...2'
cd8dev_counts_exon = cd8dev_counts_exon[, c("ACD81", "NCD81", 
                                            "ACD82", "NCD82", 
                                            "ACD83", "NCD83", 
                                            "ACD84", "NCD84")] 
cd8dev_counts_exon$gene = rownames(cd8dev_counts_exon)


################## datasets loaded ##################

# combining data sets
counts_exon_all = Reduce(function(x,y) merge(x = x, y = y, by = 'gene'), 
                         list(mir29_counts_exon, erin2015_counts_exon, blood2016_counts_exon, 
                              cleanCd8_counts_exon, mir29cre_counts_exon, cd8dev_counts_exon)) 
counts_exon_all_noname = counts_exon_all
rownames(counts_exon_all_noname) = counts_exon_all_noname$gene
counts_exon_all_noname$gene = NULL


sample_info = data.frame(colnames(counts_exon_all_noname), 
                         c(rep('Mir29 Project', 12),
                           rep('Wissink 2015', 2),
                           rep('Wang 2016', 4-2), 
                           rep('Clean Project', 12), 
                           rep('Mir29Cre Project', 12-4), 
                           rep('CD8Dev Project', 16-8)),
                         c('Adult', 'Adult', 'Neo', 'Neo', 'Adult', 'Adult', 'Neo', 'Neo', 'Adult', 'Adult', 'Neo', 'Neo', 
                           rep(c('Adult', 'Neo'), 1+2-1),
                           rep('Adult', 12), 
                           rep('Adult', 12-4), # 7+
                           "Adult", "Neo", "Adult", "Neo", "Adult", "Neo", "Adult", "Neo"),
                         c(rep(c('TN', 'MP'), 6), 
                           rep('Bulk', 6-2),
                           'TN', 'TN', 'TN', 'MP', 'MP', 'MP', 'TN', 'TN', 'TN', 'MP', 'MP', 'MP',
                           #rep('Bulk', 7),
                           rep(c('TN', 'MP'), 6-2),  
                           rep('Bulk', 16-8)),
                         c(rep('Spleen', 12), rep('Spleen', 2+2), rep('Spleen', 12), rep('Spleen', 12-4), 
                           rep(c("Spleen"), 8)), 
                         c(rep('WT', 12+2+4-2), rep('WT', 12), # +7
                           'WT', 'WT', 'Creneg', 'Creneg', 'WT', 'WT', 'Creneg', 'Creneg', rep('WT', 16-8)), 
                         c(rep('Clean', 12+2+4-2+6), rep('Dirty', 6), rep('Clean', 12-4+16-8)) 
)
colnames(sample_info) = c("Sample", "Source", "Age", "CellType", "Organ", "Treatment", "Clean") # 
write.table(sample_info, paste0(outputdir, "sample_info.txt"), quote = F, sep = "\t")


# get sample_info
sample_info = read.table('inputs/exonIntron_matrix/sample_info.txt',
                         sep = '\t')
table(sample_info$Sample == colnames(counts_exon_all_noname))

rownames(sample_info) = colnames(counts_exon_all_noname)

# get vst
sample_info$Source <- factor(sample_info$Source) 
sample_info$Age <- factor(sample_info$Age) 
sample_info$CellType <- factor(sample_info$CellType) 
sample_info$Organ <- factor(sample_info$Organ) 
sample_info$Treatment <- factor(sample_info$Treatment) 
sample_info$Clean <- factor(sample_info$Clean) 

# Define experimental factors and design matrix. 
batch <- factor(sample_info$Source) 
age <- factor(sample_info$Age) 
cellType <- factor(sample_info$CellType) 
organ <- factor(sample_info$Organ) 
treatment <- factor(sample_info$Treatment) 
clean <- factor(sample_info$Clean) 

design <- model.matrix(~age+cellType+treatment+clean) 

library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = counts_exon_all_noname,
                              colData = sample_info,
                              design= ~ Age + Treatment + Clean + CellType) # Source +   (or not full rank)
keep <- rowSums(counts(dds) > 0) >= 1
dds <- dds[keep,]


# DESeq2 vst blind=T
vsd <- vst(dds, blind=T)
counts = assay(vsd)

library(ggfortify)
ntop <- 500
library(matrixStats)
rv <- rowVars(as.matrix(counts))
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
mat <- t( counts[select, ] )

pca<-prcomp(mat)
rownames(pca$x) =colnames(counts)
autoplot(pca, label = T, label.size = 3.5) 
# PC1: 35.13%, PC2: 24.74%

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
  xlab("PC1 (35.13%)") +
  ylab("PC2 (24.74%)") +
  labs(color='Source') +
  labs(shape='Age') +
  guides(color = guide_legend(order = 1),
         shape = guide_legend(order = 2)) +
  scale_color_manual(limits = c('Clean Project', 'Mir29 Project', 'Wang 2016', 'Wissink 2015', 
                                'Mir29Cre Project', 'CD8Dev Project'), 
                     labels = c('Tabilas et al. (2022)', 'Yee Mon et al. (2021)', 'Wang et al. (2016)', 'Wissink et al. (2015)', 
                                'Mir29Cre Project', 'CD8Dev Project'), 
                     values = mycolors) #+ scale_shape_discrete(labels = c('Bulk', 'VM', 'TN'))
ggsave("RNAseq_PCA/PCA_vst_blindT_sourceAge.pdf", plot = last_plot(), device = "pdf",
       width = 6.5, height = 4, dpi = 300)



# combat
design <- model.matrix(~age+treatment+clean) # +cellType
library(sva)
counts_combat = ComBat(dat=as.matrix(counts), batch=batch, mod=design, par.prior=TRUE, prior.plots=FALSE)

library(ggfortify)
ntop <- 500
library(matrixStats)
rv <- rowVars(as.matrix(counts_combat))
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
mat <- t( counts_combat[select, ] )

pca<-prcomp(mat)
rownames(pca$x) =colnames(counts_combat)
autoplot(pca, label = T, label.size = 3.5) 
# PC1: 56.45%, PC2: 15.69%
autoplot(pca, y = 3, label = T, label.size = 3.5) 
# PC3: 7.58%

pca_plot = as.data.frame(pca$x)
pca_plot$Sample = sample_info$Sample
pca_plot$Source = sample_info$Source
pca_plot$Age = sample_info$Age
pca_plot$CellType = sample_info$CellType
pca_plot$Organ = sample_info$Organ
pca_plot$Treatment = sample_info$Treatment
pca_plot$Clean = sample_info$Clean

ggplot(pca_plot, aes(x = PC1, y= PC2)) +
  geom_point(size = 3, aes(colour = CellType, shape = Age)) +
  theme_bw(base_size = 14) +
  xlab("PC1 (56.45%)") +
  ylab("PC2 (15.69%)") +
  labs(color='Cell type') +
  labs(shape='Age') +
  guides(color = guide_legend(order = 1),
         shape = guide_legend(order = 2)) +
  #theme(legend.box = "horizontal") +
  scale_color_manual(values = mycolors[5:7], 
                     limits = c('TN', 'Bulk', 'MP'),
                     labels = c('TN', 'Bulk', 'VM'))
ggsave("RNAseq_PCA/PCA_vst_blindT_combat_celltypeAge.pdf", plot = last_plot(), device = "pdf",
       width = 5.5, height = 4, dpi = 300)


library(PCAtools)
# repeat previous result
p0 <- pca(t(mat))
biplot(p0)
# add metadata
metadata = sample_info
metadata$Sample = NULL
metadata$Source = factor(metadata$Source, levels = unique(metadata$Source)) # randomly ordered
metadata$Age = factor(metadata$Age, levels = c("Adult", 'Neo')) 
metadata$CellType = factor(metadata$CellType, levels = c("TN", 'Bulk','MP')) 
#metadata$Organ = factor(metadata$Organ, levels = c("Spleen", 'Thymus')) 
metadata$Treatment = factor(metadata$Treatment, levels = c("WT", 'CreNeg', 'iLin28b','CrePos')) 
metadata$Clean = factor(metadata$Clean, levels = c("Clean", 'Dirty')) 

metadata_tmp = metadata
colnames(metadata_tmp) = c('Source', 'Neonatal/Adult', 'VM/TN', 'Dirty/Clean') #'Thymus/Spleen', 
metadata_tmp = metadata_tmp[,rev(c(3,2))]
p <- pca(t(mat), metadata = metadata_tmp)
pdf("RNAseq_PCA/PCA_eigencorplot.pdf", width = 7, height = 4)
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
