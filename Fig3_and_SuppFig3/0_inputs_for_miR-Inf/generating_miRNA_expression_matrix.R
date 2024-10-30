# generating the following files:
# miRNA expression matrix: counts_vst_blindT_combat_expressed1000.txt
# expressed miRNAs: mirnas_expressed1000.txt

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



# cleanCd8 after removing outlier (wo: without outlier)
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

colnames(combined_raw) = c("adult_vm_1", "adult_vm_2", "adult_tn_1", "adult_tn_2", "neo_vm_1", "neo_vm_2", "neo_tn_1",
                           "adult_naive_ilm", "adult_naive_neb", "neo_naive_ilm", "neo_naive_neb", "splenic_adult", "splenic_neonate", 
                           "Clean_TN_rep2", "Clean_TN_rep3",  "Clean_VM_rep1", "Clean_VM_rep2", "Clean_VM_rep3","Dirty_TN_rep1", "Dirty_TN_rep2", "Dirty_TN_rep3", "Dirty_VM_rep1", "Dirty_VM_rep2", "Dirty_VM_rep3", # "Clean_TN_rep1", 
                           'mm1', 'mm2', 'mt1', 'mt2', 'wm1', 'wm2', 'wt1', 'wt2', # 'pm1', 'pm2', 'pt1', 'pt2', 
                           "ACD81", "NCD81", "ACD82", "NCD82", "ACD83", "NCD83", "ACD84", "NCD84"
) #"newborn_naive_ilm", "splenic_newborn", "thymic_newborn", 

write.table(combined_raw, paste0(outputdir, "counts_raw_all.txt"), quote = F, sep = "\t")


colnames(combined_norm) = c("adult_vm_1", "adult_vm_2", "adult_tn_1", "adult_tn_2", "neo_vm_1", "neo_vm_2", "neo_tn_1",
                           "adult_naive_ilm", "adult_naive_neb", "neo_naive_ilm", "neo_naive_neb", "splenic_adult", "splenic_neonate", 
                           "Clean_TN_rep2", "Clean_TN_rep3",  "Clean_VM_rep1", "Clean_VM_rep2", "Clean_VM_rep3","Dirty_TN_rep1", "Dirty_TN_rep2", "Dirty_TN_rep3", "Dirty_VM_rep1", "Dirty_VM_rep2", "Dirty_VM_rep3", # "Clean_TN_rep1", 
                           'mm1', 'mm2', 'mt1', 'mt2', 'wm1', 'wm2', 'wt1', 'wt2', # 'pm1', 'pm2', 'pt1', 'pt2', 
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


# sample_info before taking mean for cleanCd8 - 63 samples (1 outlier in cleanThy - ct1 and 1 outlier in cleanCd8 - ds3)
sample_info = data.frame(colnames(combined_raw), 
                         c(rep('Wang 2016 (unpublished)', 7), rep('Wissink 2015', 8-2), rep('Clean CD8', 11),  rep('Mir29Cre', 12-4), rep('CD8 Dev', 16-8)), # rep('Wang 2016', 2), rep('Clean Thy', 7),
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
                           0.264, 0.225, 0.49, 0.329, 0.333, 0.284, 0.39, 0.373), # 2+ mean of replicates (not including Clean_TN_rep1), # 0.593, # 0.342, 0.222, 0.418, 0.391, 0.285, 0.355, 0.421, 0.314, 
                         c(rep('Clean', 7+8-2+5), rep('Dirty', 6), rep('Clean', 12-4+16-8)) # rep('Clean', 4), rep('Dirty', 3), 
                        )
colnames(sample_info) = c("Sample", "Source", "Age", "CellType", "Organ", "Treatment", "Kit", "MappingRate", "Clean") # 
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
                              design= ~ Source + Age + Kit + Clean + Treatment ) #   + CellType

vsd = varianceStabilizingTransformation(dds, blind = TRUE)

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

# for cleanCd8 project, smRNAseq and RNAseq is not matched. taking mean across replicates
ct_countscombat = rowMeans(countscombat[,c("Clean_TN_rep2", "Clean_TN_rep3")])
cv_countscombat = rowMeans(countscombat[,c("Clean_VM_rep1", "Clean_VM_rep2", "Clean_VM_rep3")])
dt_countscombat = rowMeans(countscombat[,c("Dirty_TN_rep1", "Dirty_TN_rep2", "Dirty_TN_rep3")])
dv_countscombat = rowMeans(countscombat[,c("Dirty_VM_rep1", "Dirty_VM_rep2", "Dirty_VM_rep3")])
drops = c("Clean_TN_rep2", "Clean_TN_rep3",  "Clean_VM_rep1", "Clean_VM_rep2", "Clean_VM_rep3","Dirty_TN_rep1", "Dirty_TN_rep2", "Dirty_TN_rep3", "Dirty_VM_rep1", "Dirty_VM_rep2", "Dirty_VM_rep3")
countscombat_mean = countscombat[ , !(colnames(countscombat) %in% drops)]
countscombat_mean$Clean_TN_rep1 = ct_countscombat
countscombat_mean$Clean_TN_rep2 = ct_countscombat
countscombat_mean$Clean_TN_rep3 = ct_countscombat
countscombat_mean$Clean_VM_rep1 = cv_countscombat
countscombat_mean$Clean_VM_rep2 = cv_countscombat
countscombat_mean$Clean_VM_rep3 = cv_countscombat
countscombat_mean$Dirty_TN_rep1 = dt_countscombat
countscombat_mean$Dirty_TN_rep2 = dt_countscombat
countscombat_mean$Dirty_TN_rep3 = dt_countscombat
countscombat_mean$Dirty_VM_rep1 = dv_countscombat
countscombat_mean$Dirty_VM_rep2 = dv_countscombat
countscombat_mean$Dirty_VM_rep3 = dv_countscombat
countscombat_mean = countscombat_mean[, c("adult_vm_1", "adult_vm_2", "adult_tn_1", "adult_tn_2", "neo_vm_1", "neo_vm_2", "neo_tn_1",
                                          "adult_naive_ilm", "adult_naive_neb", "neo_naive_ilm", "neo_naive_neb", "splenic_adult", "splenic_neonate", 
                                          "Clean_TN_rep1", "Clean_TN_rep2", "Clean_TN_rep3",  "Clean_VM_rep1", "Clean_VM_rep2", "Clean_VM_rep3","Dirty_TN_rep1", "Dirty_TN_rep2", "Dirty_TN_rep3", "Dirty_VM_rep1", "Dirty_VM_rep2", "Dirty_VM_rep3",
                                          'mm1', 'mm2', 'mt1', 'mt2', 'wm1', 'wm2', 'wt1', 'wt2', 
                                          "ACD81", "NCD81", "ACD82", "NCD82", "ACD83", "NCD83", "ACD84", "NCD84")]
write.table(countscombat_mean, paste0(outputdir, "counts_vst_blindT_combat.txt"), quote = F, sep = "\t")


counts_combat_expressed1000 = countscombat_mean[expressedmirnas1000, ]
write.table(counts_combat_expressed1000, paste0(outputdir, "counts_vst_blindT_combat_expressed1000.txt"), quote = F, sep = "\t")
