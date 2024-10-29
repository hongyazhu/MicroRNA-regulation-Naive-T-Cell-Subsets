# Preprocessing, batch effect removal and PCA-related analysis for small RNAseq datasets.
# Fig 1A-C, Supplemental Fig 1B.

library(edgeR) 
library(ggfortify)
library(limma) 
library(DESeq2)
library(ggplot2)
library(ggrepel)

outputdir = 'inputs/mirna_exp_matrix/'


#########################################################################
### read data in (all mirna families), raw and norm counts seperately ###
#########################################################################

blood2016unpub = read.csv("data_analysis/adultNeoTNVM-smRNAseq/all_family_all.txt", header = T, sep = "\t")
blood2016unpub$total = NULL
blood2016unpubsamplesize = (ncol(blood2016unpub) - 1)/2
blood2016unpub_raw = blood2016unpub[, c(1:(1+blood2016unpubsamplesize))]
blood2016unpub_norm = blood2016unpub[, c(1, (2+blood2016unpubsamplesize):ncol(blood2016unpub))]

erin2015 = read.csv("data_analysis/Wissink2015_smRNAseq_all/all_family_all.txt", header = T, sep = "\t")
erin2015$total = NULL
erin2015_nonbn_noeff = erin2015
erin2015_nonbn_noeff$nbn = NULL
erin2015_nonbn_noeff$nbn.norm. = NULL
erin2015_nonbn_noeff$snb = NULL
erin2015_nonbn_noeff$snb.norm. = NULL
erin2015_nonbn_noeff$tnb = NULL
erin2015_nonbn_noeff$tnb.norm. = NULL

erin2015_nonbn_noeff$a5i = NULL
erin2015_nonbn_noeff$a5i.norm. = NULL
erin2015_nonbn_noeff$a5n = NULL
erin2015_nonbn_noeff$a5n.norm. = NULL
erin2015_nonbn_noeff$a7i = NULL
erin2015_nonbn_noeff$a7i.norm. = NULL
erin2015_nonbn_noeff$a7n = NULL
erin2015_nonbn_noeff$a7n.norm. = NULL

erin2015_nonbn_noeff$n5i = NULL
erin2015_nonbn_noeff$n5i.norm. = NULL
erin2015_nonbn_noeff$n5n = NULL
erin2015_nonbn_noeff$n5n.norm. = NULL
erin2015_nonbn_noeff$n7i = NULL
erin2015_nonbn_noeff$n7i.norm. = NULL
erin2015_nonbn_noeff$n7n = NULL
erin2015_nonbn_noeff$n7n.norm. = NULL

erin2015_nonbn_noeff$tad = NULL
erin2015_nonbn_noeff$tad.norm. = NULL
erin2015_nonbn_noeff$tnn = NULL
erin2015_nonbn_noeff$tnn.norm. = NULL

erin2015_nonbn_noeffsamplesize = (ncol(erin2015_nonbn_noeff) - 1)/2
erin2015_nonbn_noeff_raw = erin2015_nonbn_noeff[, c(1:(1+erin2015_nonbn_noeffsamplesize))]
erin2015_nonbn_noeff_norm = erin2015_nonbn_noeff[, c(1, (2+erin2015_nonbn_noeffsamplesize):ncol(erin2015_nonbn_noeff))]


# cleanCd8, removing outlier (wo: without outlier)
cleanCd8wo = read.csv("data_analysis/cleandirtyCD8_smRNAseq/May202021/all_family_all.txt", header = T, sep = "\t")
cleanCd8wo$total = NULL
# remove outlier
cleanCd8wo$ct1 = NULL
cleanCd8wo$ct1.norm. = NULL
cleanCd8wosamplesize = (ncol(cleanCd8wo) - 1)/2
cleanCd8wo_raw = cleanCd8wo[, c(1:(1+cleanCd8wosamplesize))]
cleanCd8wo_norm = cleanCd8wo[, c(1, (2+cleanCd8wosamplesize):ncol(cleanCd8wo))]


mir29cre = read.csv("data_analysis/mir29cre-smRNAseq-full/all_family_all.txt", header = T, sep = "\t")
mir29cre$total = NULL
mir29cre$pm1 = NULL
mir29cre$pm1.norm. = NULL
mir29cre$pm2 = NULL
mir29cre$pm2.norm. = NULL
mir29cre$pt1 = NULL
mir29cre$pt1.norm. = NULL
mir29cre$pt2 = NULL
mir29cre$pt2.norm. = NULL
mir29cresamplesize = (ncol(mir29cre) - 1)/2
mir29cre_raw = mir29cre[, c(1:(1+mir29cresamplesize))]
mir29cre_norm = mir29cre[, c(1, (2+mir29cresamplesize):ncol(mir29cre))]

cd8dev = read.csv("data_analysis/cd8dev_smallRNAseq/all_family_all.txt", header = T, sep = "\t")
cd8dev$total = NULL
cd8devsamplesize = (ncol(cd8dev) - 1)/2
cd8dev_raw = cd8dev[, c(1:(1+cd8devsamplesize))]
cd8dev_norm = cd8dev[, c(1, (2+cd8devsamplesize):ncol(cd8dev))]
# get cd8dev names:
library("readxl")
filenames <- read_excel('filenames.xlsx', col_names = F)
filenames$num = paste0('00', filenames$...1)
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}
filenames$num = paste0('X', substrRight(filenames$num, 3))
colnames(cd8dev_raw) = as.data.frame(filenames[match(colnames(cd8dev_raw),filenames$num),'...2'])$'...2'
colnames(cd8dev_raw)[1] = 'seed'
colnames(cd8dev_norm) = paste0(as.data.frame(filenames[match(substr(colnames(cd8dev_norm), 1, 4),filenames$num),'...2'])$'...2', '.norm.')
colnames(cd8dev_norm)[1] = 'seed'
cd8dev_raw_laterstages = cd8dev_raw[ , grepl( "seed|CD8" , names( cd8dev_raw ) ) ] # SP|
cd8dev_norm_laterstages = cd8dev_norm[ , grepl( "seed|CD8" , names( cd8dev_norm ) ) ] # SP|


seed_mirna_all = data.frame(miRNA = rownames(blood2016unpub),
                            seed = blood2016unpub$seed)


###############
### combine ###
###############

# combine data frames - raw
tmp2 = merge(x = erin2015_nonbn_noeff_raw, y = cleanCd8wo_raw, by = 'seed', all = TRUE)
tmp3 = merge(x = blood2016unpub_raw, y = tmp2, by = "seed", all = TRUE)
tmp5 = merge(x = tmp3, y = mir29cre_raw, by = "seed", all = TRUE)
combined_raw = merge(x = tmp5, y = cd8dev_raw_laterstages, by = "seed", all = TRUE)
row.names(combined_raw) <- combined_raw$seed
combined_raw$seed <- NULL
# select rows with at least 1 positive vluaes
combined_raw = combined_raw[rowSums(combined_raw)>0,]

# combine data frames - norm
tmp2 = merge(x = erin2015_nonbn_noeff_norm, y = cleanCd8wo_norm, by = 'seed', all = TRUE)
tmp3 = merge(x = blood2016unpub_norm, y = tmp2, by = "seed", all = TRUE)
tmp5 = merge(x = tmp3, y = mir29cre_norm, by = "seed", all = TRUE)
combined_norm = merge(x = tmp5, y = cd8dev_norm_laterstages, by = "seed", all = TRUE)
row.names(combined_norm) <- combined_norm$seed
combined_norm$seed <- NULL
# select rows with at least 1 positive vluaes
combined_norm = combined_norm[rowSums(combined_norm)>0,]

colnames(combined_raw) = c(
  "adult_vm_1", "adult_vm_2", "adult_tn_1", "adult_tn_2", "neo_vm_1", "neo_vm_2", "neo_tn_1",
  "adult_naive_ilm", "adult_naive_neb", "neo_naive_ilm", "neo_naive_neb", "splenic_adult", "splenic_neonate", 
  "Clean_TN_rep2", "Clean_TN_rep3",  "Clean_VM_rep1", "Clean_VM_rep2", "Clean_VM_rep3","Dirty_TN_rep1", "Dirty_TN_rep2", "Dirty_TN_rep3", "Dirty_VM_rep1", "Dirty_VM_rep2", "Dirty_VM_rep3", # "Clean_TN_rep1", 
  'mm1', 'mm2', 'mt1', 'mt2', 'wm1', 'wm2', 'wt1', 'wt2', 
  "ACD81", "NCD81", "ACD82", "NCD82", "ACD83", "NCD83", "ACD84", "NCD84"
) 

write.table(combined_raw, paste0(outputdir, "counts_raw_all.txt"), quote = F, sep = "\t")
write.table(rownames(combined_raw), "../mirna_DE/all_mirnas_nonbn.txt", row.names = FALSE, col.names=FALSE, sep = "\t", quote = F)


colnames(combined_norm) = c(
  "adult_vm_1", "adult_vm_2", "adult_tn_1", "adult_tn_2", "neo_vm_1", "neo_vm_2", "neo_tn_1",
  "adult_naive_ilm", "adult_naive_neb", "neo_naive_ilm", "neo_naive_neb", "splenic_adult", "splenic_neonate", 
  "Clean_TN_rep2", "Clean_TN_rep3",  "Clean_VM_rep1", "Clean_VM_rep2", "Clean_VM_rep3","Dirty_TN_rep1", "Dirty_TN_rep2", "Dirty_TN_rep3", "Dirty_VM_rep1", "Dirty_VM_rep2", "Dirty_VM_rep3", # "Clean_TN_rep1", 
  'mm1', 'mm2', 'mt1', 'mt2', 'wm1', 'wm2', 'wt1', 'wt2', 
  "ACD81", "NCD81", "ACD82", "NCD82", "ACD83", "NCD83", "ACD84", "NCD84"
) 

write.table(combined_norm, paste0(outputdir, "combined_norm_all.txt"), quote = F, sep = "\t")


# filter mirnas lowly expressed in all datasets using combined_norm - input for batch effect removal methods
combined_norm_expressed10 = combined_norm[rowSums(combined_norm > 10) >= 1, ]
expressedmirnas10 = rownames(combined_norm_expressed10)
combined_raw_expressed10 = combined_raw[expressedmirnas10, ]


# keep mirnas highly expressed in at least one sample using combined_norm
combined_norm_expressed1000 = combined_norm[rowSums(combined_norm > 1000) >= 1, ]
expressedmirnas1000 = rownames(combined_norm_expressed1000)
write.table(expressedmirnas1000, paste0(outputdir, "../mirna_DE/mirnas_expressed1000.txt"), row.names = FALSE, col.names=FALSE, sep = "\t", quote = F)
combined_raw_expressed1000 = combined_raw[expressedmirnas1000, ]
write.table(combined_raw_expressed1000, paste0(outputdir, "counts_raw_expressed1000.txt"), quote = F, sep = "\t")

pca <- prcomp(t(combined_norm_expressed1000))
autoplot(pca, label = T, label.size = 3.5, label.repel = T) 

sample_info = data.frame(colnames(combined_raw), 
                         c(rep('Wang 2016 (unpublished)', 7), rep('Wissink 2015', 8-2), rep('Clean CD8', 11),  rep('Mir29Cre', 12-4), rep('CD8 Dev', 16-8)), 
                         c(#'Adult', 'Adult', 
                           rep('Adult', 4), rep('Neo', 3),
                           rep('Adult', 2), rep('Neo', 2), rep(c('Adult', 'Neo'), 2-1),
                           rep('Adult', 11), 
                           #rep('Adult', 7), 
                           rep('Adult', 12-4),
                           rep(c('Adult', 'Neo', 'Adult', 'Neo'), 4-2)),
                         c(#'All', 'All', 
                           'MP', 'MP', 'TN', 'TN', 'MP', 'MP', 'TN', 
                           rep('All', 8-2),
                           'TN', 'TN', 'MP', 'MP', 'MP', 'TN', 'TN', 'TN', 'MP', 'MP', 'MP', # 'TN', 
                           #rep('All', 7),
                           'MP', 'MP', 'TN', 'TN', 'MP', 'MP', 'TN', 'TN', # 'MP', 'MP', 'TN', 'TN', 
                           rep('All', 16-8)),
                         c(#rep('Spleen', 2), 
                           rep('Spleen', 7), rep('Spleen', 2+2+2), rep('Spleen', 11), rep('Spleen', 12-4), rep(c('Spleen'), 8)), 
                         c(#'ilin28b', 1+
                           rep('WT', 7+8-2), rep('WT', 11), 'Creneg', 'Creneg', 'Creneg', 'Creneg', 'WT', 'WT', 'WT', 'WT', rep('WT', 16-8)), 
                         c(#rep('NEB', 2),
                           rep('Illumina', 7),
                           rep(c('Illumina', 'NEB'), 2), rep('Illumina', 4-2),
                           rep('NEB', 11+12-4+16-8)), #+7
                         c(#0.122, 0.071,
                           0.475, 0.545, 0.526, 0.511, 0.446, 0.446, 0.469,
                           0.647, 0.417, 0.481, 0.226, 0.629, 0.461, # 0.692, 0.637, 
                           0.642, 0.628, 0.686, 0.69, 0.684, 0.744, 0.707, 0.66, 0.655, 0.593, 0.638, # 0.599, 
                           #0.608, 0.707, 0.67, 0.714, 0.594, 0.642, 0.66,
                           0.376, 0.457, 0.559, 0.605, 0.477, 0.542, 0.548, 0.457, # , 0.323, 0.612, 0.495, 0.599
                           0.264, 0.225, 0.49, 0.329, 0.333, 0.284, 0.39, 0.373), 
                         c(rep('Clean', 7+8-2+5), rep('Dirty', 6), rep('Clean', 12-4+16-8))
)
colnames(sample_info) = c("Sample", "Source", "Age", "CellType", "Organ", "Treatment", "Kit", "MappingRate", "Clean") 
write.table(sample_info, paste0(outputdir, "sample_info.txt"), quote = F, sep = "\t")


###########
### vst ###
###########

# vst_blindT, using combined_raw_expressed10 as input (count matrix filtered mirnas not lowly expressed in all samples)

# DESeqDataSetFromMatrix does not allow non-integer values, so convert the matrix to integer first
combined_int = combined_raw_expressed10
for (i in 1:nrow(combined_raw_expressed10))
  for(j in 1:ncol(combined_raw_expressed10))
    combined_int[i,j] = as.integer(combined_raw_expressed10[i,j])
dds <- DESeqDataSetFromMatrix(countData = combined_int, 
                              colData = sample_info,
                              design= ~ Source + Age  + Kit + Clean + Treatment ) #   ++ CellType 

# DESeq2 vst blind=T
vsd = varianceStabilizingTransformation(dds, blind = TRUE)

### pca for vst

pcaData <- plotPCA(vsd, intgroup=c("Source", "Age", "Organ", "Treatment", "Kit", "Clean"), returnData=TRUE) #"CellType", , "State"
percentVar <- round(100 * attr(pcaData, "percentVar"))
library(RColorBrewer)
mycolors = c(brewer.pal(name="Set2", n = 7))
ggplot(pcaData, aes(PC1, PC2, color=Source, shape=Kit)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_bw(base_size = 14, base_family='ArialMT') +
  scale_color_manual(limits = c("Wissink 2015", "Wang 2016 (unpublished)", "Mir29Cre", "Clean CD8", "CD8 Dev"), # , "Clean Thy"
                     labels = c('Wissink et al. (2015)', 'Unpublished (Wang 2016)', 'Unpublished (Mir29Cre)', 'Unpublished (Clean CD8)', 'Unpublished (CD8Dev)'), # 'Unpublished (CleanThy)', 
                     values = mycolors)+
  guides(color = guide_legend(order = 1),
         shape = guide_legend(order = 2))
ggsave(
  paste0("smRNAseq_PCA/PCA_vst_blindT_sourceKit.pdf"),
  plot = last_plot(),
  device = "pdf",
  width = 6.5,
  height = 4,
  dpi = 300
)

assayvsd = as.data.frame(assay(vsd))


##############
### combat ###
##############

# combat 
library(sva)

# Define experimental factors and design matrix. 
batch <- factor(sample_info$Source) 
age <- factor(sample_info$Age) 
organ <- factor(sample_info$Organ) 
treatment <- factor(sample_info$Treatment) 
kit <- factor(sample_info$Kit) 
clean <- factor(sample_info$Clean) 
celltype <- factor(sample_info$CellType) 


batch_kit <- paste(kit, batch, sep = "_")
design <- model.matrix(~age+clean+treatment) 
countscombat = ComBat(dat=as.matrix(assayvsd), batch=batch_kit, mod=design, par.prior=TRUE, prior.plots=FALSE)
countscombat = as.data.frame(countscombat)

### pca for combat

library(ggfortify)
ntop <- 500
rv <- rowVars(countscombat)
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
mat <- t( countscombat[select, ] )

pca<-prcomp(mat)
rownames(pca$x) =colnames(countscombat)
autoplot(pca, label = T, label.size = 3.5) 
autoplot(pca, label = T, label.size = 3.5, y = 3) 

pca_plot = as.data.frame(pca$x)
pca_plot$Sample = sample_info$Sample
pca_plot$Source = sample_info$Source
pca_plot$Age = sample_info$Age
pca_plot$CellType = sample_info$CellType
pca_plot$Organ = sample_info$Organ
pca_plot$Treatment = sample_info$Treatment
pca_plot$Kit = sample_info$Kit
pca_plot$MappingRate = sample_info$MappingRate
pca_plot$Clean = sample_info$Clean


assay(vsd) <- countscombat ### !! assay(vsd) is changed to combat value
pcaData <- plotPCA(vsd, intgroup=c("Source", "Age", "Organ", "Treatment", "Kit", "CellType", "MappingRate", "Clean"), returnData=TRUE) # "State", 
percentVar <- round(100 * attr(pcaData, "percentVar"), digits = 2)

ggplot(pcaData, aes(PC1, PC2, color=age, shape=CellType)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_bw(base_size = 14, base_family='ArialMT')+
  scale_color_manual(values = mycolors[3:4])+
  guides(color = guide_legend(order = 1, title = 'Age'),
         shape = guide_legend(order = 2, title = 'CellType'))
ggsave(
  paste0("smRNAseq_PCA/PCA_vst_blindT_combat_ageCellType.pdf"),
  plot = last_plot(),
  device = "pdf",
  width =5.5,
  height = 4,
  dpi = 300
)


library(PCAtools)
# repeat previous result
# add metadata
metadata = sample_info
rownames(metadata) = gsub('-', '.', metadata$Sample)
metadata$Sample = NULL
metadata$Source = factor(metadata$Source, levels = unique(metadata$Source)) # randomly ordered
metadata$Age = factor(metadata$Age, levels = c("Adult", 'Neo')) 
metadata$CellType = factor(metadata$CellType, levels = c("TN", 'All','MP')) 
metadata$Clean = factor(metadata$Clean, levels = c("Clean", 'Dirty')) 

metadata_tmp = metadata
colnames(metadata_tmp) = c('Source', 'Neonatal/Adult', 'VM/TN', 'Organ', 'Dirty/Clean') # 'Thymus/Spleen', 
metadata_tmp = metadata_tmp[,rev(c(3,2))]
p <- pca(countscombat, metadata = metadata_tmp)
biplot(p)
pdf("smRNAseq_PCA/PCA_combat_eigencorplot.pdf", width = 7, height = 4)
eigencorplot(p,
             components = getComponents(p, 1:3),
             metavars = rev(c('VM/TN','Neonatal/Adult')), # ,'Thymus/Spleen'
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



# get loadings
autoplot(pca, label = T, label.size = 3.5) 
loadings <- as.data.frame(pca$rotation)[,1:3]

# get miRNA names 
mir_family_info <- read.table("/path/to/targetscan/mouse/miR_Family_Info.txt", head = T, fill = T, sep = "\t") 
hsa_mir_family_info <- mir_family_info[substr(mir_family_info[,'MiRBase.ID'] ,1 ,3) == "mmu", ] 
hsa_con_mir_family_info <- hsa_mir_family_info[hsa_mir_family_info[,'Family.Conservation.'] >= 2, ] 
seed_con <- hsa_con_mir_family_info[,c('Seed.m8', 'MiRBase.ID')] 
mirna_seeds_con <- aggregate(seed_con[,2], by = list(seed = seed_con$Seed), FUN = paste) 
colnames(mirna_seeds_con) <- c("seed", "miRNA") 
mirna_seeds_con$miRNA_name = c('miR-106ab,17,20ab,93-5p,6383','miR-141,200a-3p','miR-132,212-3p','miR-451a','miR-490-3p','miR-191-5p','miR-124-3p.1','miR-18ab-5p','miR-291a,294,295,302abd-3p','miR-200bc,429-3p','miR-216a-5p','miR-216b-5p','miR-365-3p','miR-101a-3p.1','miR-144-3p','miR-181abcd-5p','miR-140-3p.2,miR-497b','miR-100,99ab-5p','miR-10ab-5p','miR-217-5p','miR-193ab-3p','miR-383-5p.2','miR-29abc-3p','miR-15ab,16,1907,195a,322,497a-5p,miR-195b,6342,6353,6419','miR-129-1,2-3p','miR-22-3p','miR-21a,590-5p,miR-21c','miR-196ab-5p','miR-130ab,301ab-3p,miR-130c,6341,6389,721','miR-302c-3p','miR-140-5p,miR-876-3p','miR-142a-5p','miR-425-5p,miR-489-3p','miR-183-5p','miR-135ab-5p','miR-455-5p,miR-5129-3p','miR-25,363,367,92ab-3p,miR-32-5p','miR-128-3p,miR-6539','miR-802-5p','miR-199ab-3p','miR-455-3p.1','miR-148ab,152-3p','miR-140-3p.1','miR-338-3p','miR-199ab-5p','miR-125ab,351-5p,miR-6367,6394','miR-205-5p','miR-212-5p','miR-551b-3p','miR-126a-3p.1','miR-187-3p','miR-139-5p','miR-150-5p,miR-5127','miR-9-5p','miR-203-3p.1','miR-146ab-5p','miR-143-3p','let-7abcdefgik,98-5p,miR-1961','miR-190ab-5p','miR-383-5p.1','miR-219a-5p','miR-103,107-3p','miR-214-5p','miR-221,222-3p,miR-1928','miR-138-5p','miR-7ab-5p','miR-1a,206-3p,miR-1957b,6349,6382','miR-184-3p','miR-122-5p','miR-31-5p','miR-34abc,449ac-5p,-miR-449b','miR-24-3p,miR-5124b,6361,6369,6410,6413','miR-193a-5p','miR-30abcde,384-5p','miR-194-5p','miR-126a-3p.2','miR-142a-3p.1','miR-223-3p','miR-19ab-3p','miR-208ab-3p','miR-499-5p','miR-124,5624-3p,miR-6540-5p','miR-155-5p','miR-101ab-3p','miR-142a-3p.2','miR-137-3p','miR-26ab-5p','miR-27ab-3p','miR-23ab-3p','miR-145a-5p,miR-145b','miR-204,211-5p,miR-7670-3p','miR-202-5p','miR-203-3p.2','miR-192,215-5p','miR-455-3p.2,miR-682','miR-153-3p','miR-33-5p','miR-183-5p.2','miR-133a-3p.1','miR-147-3p','miR-210-3p','miR-218-5p,miR-7002-3p','miR-182-5p','miR-96-5p','miR-133ab-3p,miR-133c','miR-375-3p','miR-129-5p')
loadings_names = merge(loadings, mirna_seeds_con, by.x = 0, by.y = 'seed', all.x = TRUE)

write.csv(apply(loadings_names,2,as.character), 'smRNAseq_PCA/PCA_combat_loadings123.txt', quote = F, sep = '\t')

for (i in 1:nrow(loadings_names)){
  if (loadings_names[i, 'PC1'] > 0){
    loadings_names[i, 'PC1_info'] = 'PC1+'
  } else if (loadings_names[i, 'PC1'] < 0){
    loadings_names[i, 'PC1_info'] = 'PC1-'
  }
}
loadings_sortPC1 <- loadings_names[order(loadings_names$PC1),]
ngenes = 5
loadings_sortPC1_toplot <- loadings_sortPC1[c(1:ngenes, (nrow(loadings_names)-ngenes+1):nrow(loadings_names)),]
ggplot(loadings_sortPC1_toplot, aes(y=reorder(miRNA_name, PC1), x=PC1)) +
  geom_point(stat='identity', aes(color=PC1_info), size = 3, shape = 17) +
  ylab('')+
  #ggtitle(paste0(largeclusterName, "_dRNA_top", ngenes)) + 
  #scale_color_discrete() +
  theme_bw() + 
  theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
        axis.ticks = element_line(colour = "grey70", size = 0.2),
        panel.grid.major = element_line(colour = "grey70", size = 0.2),
        panel.grid.minor = element_blank())+
  geom_vline(aes(xintercept = 0)) + 
  scale_color_manual(limits = c('PC1+', 'PC1-'), values = c('black', 'black'))
ggsave(filename = paste0("smRNAseq_PCA/PC1_combat_loadings_top", ngenes, ".pdf"), 
       device = "pdf", height = 3, width = 5.5)

