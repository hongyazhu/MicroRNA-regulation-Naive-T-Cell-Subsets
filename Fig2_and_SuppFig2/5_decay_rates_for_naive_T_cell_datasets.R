# Calculating decay rates. Output:
# decay_rates.txt 
# decay_rates_match_noNA_noInf_noall0_10exon.txt (matched to smallRNAseq to be used as input of miR-Inf, filtered)


library(ggplot2)
outputdir = 'inputs/exonIntron_matrix/'

tpm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(rate) * 1e6
}




#############
### mir29 ###
#############


mir29_counts_gb = read.table("featureCounts_output/mir29_counts_gb.txt", header = T, row.names = 1)
gb_length = mir29_counts_gb$Length
mir29_counts_gb = mir29_counts_gb[,6:ncol(mir29_counts_gb)]
colnames(mir29_counts_gb) = c("NeoM29-3-VM", 
                              "Adult-1-TN", "Adult-1-VM", "NeoNC-1-TN", "NeoNC-1-VM", "NeoM29-1-TN", "NeoM29-1-VM", 
                              "Adult-2-TN", "Adult-2-VM", "NeoNC-2-TN", "NeoNC-2-VM", "NeoM29-2-TN", "NeoM29-2-VM",
                              "Adult-3-TN", "Adult-3-VM", "NeoNC-3-TN", "NeoNC-3-VM", "NeoM29-3-TN") 
mir29_counts_gb = mir29_counts_gb[, c("Adult-1-TN", "Adult-1-VM", "NeoNC-1-TN", "NeoNC-1-VM", 
                                      "Adult-2-TN", "Adult-2-VM", "NeoNC-2-TN", "NeoNC-2-VM", 
                                      "Adult-3-TN", "Adult-3-VM", "NeoNC-3-TN", "NeoNC-3-VM") ]

mir29_counts_exon = read.table("featureCounts_output/mir29_counts_exon.txt", header = T, row.names = 1)
exon_length = mir29_counts_exon$Length
mir29_counts_exon = mir29_counts_exon[,6:ncol(mir29_counts_exon)]
colnames(mir29_counts_exon) = c("NeoM29-3-VM", 
                                "Adult-1-TN", "Adult-1-VM", "NeoNC-1-TN", "NeoNC-1-VM", "NeoM29-1-TN", "NeoM29-1-VM", 
                                "Adult-2-TN", "Adult-2-VM", "NeoNC-2-TN", "NeoNC-2-VM", "NeoM29-2-TN", "NeoM29-2-VM",
                                "Adult-3-TN", "Adult-3-VM", "NeoNC-3-TN", "NeoNC-3-VM", "NeoM29-3-TN") 
mir29_counts_exon = mir29_counts_exon[, c("Adult-1-TN", "Adult-1-VM", "NeoNC-1-TN", "NeoNC-1-VM", 
                                          "Adult-2-TN", "Adult-2-VM", "NeoNC-2-TN", "NeoNC-2-VM", 
                                          "Adult-3-TN", "Adult-3-VM", "NeoNC-3-TN", "NeoNC-3-VM") ]

#table(rownames(mir29_counts_gb) == rownames(mir29_counts_exon))

mir29_counts_intron <- mir29_counts_gb - mir29_counts_exon
intron_length = gb_length - exon_length
# <0 intron counts: mir29_counts_intron[rowSums(mir29_counts_intron < 0) > 0,]
#write.table(mir29_counts_intron[rowSums(mir29_counts_intron < 0) > 0,], 'test.txt', quote = F)'

mir29_counts_intron_noNeg <- mir29_counts_intron[rowSums(mir29_counts_intron < 0) == 0,]
mir29_counts_exon_noNeg <- mir29_counts_exon[rowSums(mir29_counts_intron < 0) == 0,]

intron_length_noNeg = intron_length[rowSums(mir29_counts_intron < 0) == 0]
exon_length_noNeg = exon_length[rowSums(mir29_counts_intron < 0) == 0]


mir29_counts_intron_wIntron = mir29_counts_intron_noNeg[intron_length_noNeg != 0, ]
mir29_counts_exon_wIntron = mir29_counts_exon_noNeg[intron_length_noNeg != 0, ]

intron_length_wIntron = intron_length_noNeg[intron_length_noNeg != 0]
exon_length_wIntron = exon_length_noNeg[intron_length_noNeg != 0]


mir29_tpm_exon_wIntron <- apply(mir29_counts_exon_wIntron, 2, function(x) tpm(x, exon_length_wIntron))
mir29_tpm_intron_wIntron <- apply(mir29_counts_intron_wIntron, 2, function(x) tpm(x, intron_length_wIntron))

mir29_tpm_ex_in = merge(mir29_tpm_exon_wIntron, mir29_tpm_intron_wIntron, by = 0, suffixes = c(".exon",".intron"))
rownames(mir29_tpm_ex_in) = mir29_tpm_ex_in$Row.names
mir29_tpm_ex_in$Row.names = NULL

mir29_tpm_ex_in$gene = rownames(mir29_tpm_ex_in)


################
### erin2015 ###
################


erin2015_counts_gb = read.table("featureCounts_output/erin2015_counts_gb.txt", header = T, row.names = 1)
gb_length = erin2015_counts_gb$Length
erin2015_counts_gb = erin2015_counts_gb[,6:ncol(erin2015_counts_gb)]
colnames(erin2015_counts_gb) = gsub("X.home.hz543.data.Erin2015_RNAseq_all.Jan_2021.combine_fastq_to_sample.hisat2.", "", colnames(erin2015_counts_gb))
colnames(erin2015_counts_gb) = gsub(".bam", "", colnames(erin2015_counts_gb))
erin2015_counts_gb = erin2015_counts_gb[, c('adult_naive', 'neonate_naive')]

erin2015_counts_exon = read.table("featureCounts_output/erin2015_counts_exon.txt", header = T, row.names = 1)
exon_length = erin2015_counts_exon$Length
erin2015_counts_exon = erin2015_counts_exon[,6:ncol(erin2015_counts_exon)]
colnames(erin2015_counts_exon) = gsub("X.home.hz543.data.Erin2015_RNAseq_all.Jan_2021.combine_fastq_to_sample.hisat2.", "", colnames(erin2015_counts_exon))
colnames(erin2015_counts_exon) = gsub(".bam", "", colnames(erin2015_counts_exon))
erin2015_counts_exon = erin2015_counts_exon[, c('adult_naive', 'neonate_naive')]


#table(rownames(erin2015_counts_gb) == rownames(erin2015_counts_exon))

erin2015_counts_intron <- erin2015_counts_gb - erin2015_counts_exon
intron_length = gb_length - exon_length

erin2015_counts_intron_noNeg <- erin2015_counts_intron[rowSums(erin2015_counts_intron < 0) == 0,]
erin2015_counts_exon_noNeg <- erin2015_counts_exon[rowSums(erin2015_counts_intron < 0) == 0,]

intron_length_noNeg = intron_length[rowSums(erin2015_counts_intron < 0) == 0]
exon_length_noNeg = exon_length[rowSums(erin2015_counts_intron < 0) == 0]


erin2015_counts_intron_wIntron = erin2015_counts_intron_noNeg[intron_length_noNeg != 0, ]
erin2015_counts_exon_wIntron = erin2015_counts_exon_noNeg[intron_length_noNeg != 0, ]

intron_length_wIntron = intron_length_noNeg[intron_length_noNeg != 0]
exon_length_wIntron = exon_length_noNeg[intron_length_noNeg != 0]


erin2015_tpm_exon_wIntron <- apply(erin2015_counts_exon_wIntron, 2, function(x) tpm(x, exon_length_wIntron))
erin2015_tpm_intron_wIntron <- apply(erin2015_counts_intron_wIntron, 2, function(x) tpm(x, intron_length_wIntron))

erin2015_tpm_ex_in = merge(erin2015_tpm_exon_wIntron, erin2015_tpm_intron_wIntron, by = 0, suffixes = c(".exon",".intron"))
rownames(erin2015_tpm_ex_in) = erin2015_tpm_ex_in$Row.names
erin2015_tpm_ex_in$Row.names = NULL

erin2015_tpm_ex_in$gene = rownames(erin2015_tpm_ex_in)





#################
### blood2016 ###
#################



blood2016_counts_gb = read.table("featureCounts_output/blood2016_counts_gb.txt", header = T, row.names = 1)
gb_length = blood2016_counts_gb$Length
blood2016_counts_gb = blood2016_counts_gb[,6:ncol(blood2016_counts_gb)]
colnames(blood2016_counts_gb) = sub(".bam", "", colnames(blood2016_counts_gb))
colnames(blood2016_counts_gb) = sub("X.home.hz543.data.2016Blood.neo.adu.April032021.hisat2.", "", colnames(blood2016_counts_gb))
blood2016_counts_gb = blood2016_counts_gb[, c('adult_spleen', 'neonatal_spleen')] # , 'adult_thymus', 'neonatal_thymus'

blood2016_counts_exon = read.table("featureCounts_output/blood2016_counts_exon.txt", header = T, row.names = 1)
exon_length = blood2016_counts_exon$Length
blood2016_counts_exon = blood2016_counts_exon[,6:ncol(blood2016_counts_exon)]
colnames(blood2016_counts_exon) = sub(".bam", "", colnames(blood2016_counts_exon))
colnames(blood2016_counts_exon) = sub("X.home.hz543.data.2016Blood.neo.adu.April032021.hisat2.", "", colnames(blood2016_counts_exon))
blood2016_counts_exon = blood2016_counts_exon[, c('adult_spleen', 'neonatal_spleen')] # , 'adult_thymus', 'neonatal_thymus'



blood2016_counts_intron <- blood2016_counts_gb - blood2016_counts_exon
intron_length = gb_length - exon_length

blood2016_counts_intron_noNeg <- blood2016_counts_intron[rowSums(blood2016_counts_intron < 0) == 0,]
blood2016_counts_exon_noNeg <- blood2016_counts_exon[rowSums(blood2016_counts_intron < 0) == 0,]

intron_length_noNeg = intron_length[rowSums(blood2016_counts_intron < 0) == 0]
exon_length_noNeg = exon_length[rowSums(blood2016_counts_intron < 0) == 0]


blood2016_counts_intron_wIntron = blood2016_counts_intron_noNeg[intron_length_noNeg != 0, ]
blood2016_counts_exon_wIntron = blood2016_counts_exon_noNeg[intron_length_noNeg != 0, ]

intron_length_wIntron = intron_length_noNeg[intron_length_noNeg != 0]
exon_length_wIntron = exon_length_noNeg[intron_length_noNeg != 0]


blood2016_tpm_exon_wIntron <- apply(blood2016_counts_exon_wIntron, 2, function(x) tpm(x, exon_length_wIntron))
blood2016_tpm_intron_wIntron <- apply(blood2016_counts_intron_wIntron, 2, function(x) tpm(x, intron_length_wIntron))

blood2016_tpm_ex_in = merge(blood2016_tpm_exon_wIntron, blood2016_tpm_intron_wIntron, by = 0, suffixes = c(".exon",".intron"))
rownames(blood2016_tpm_ex_in) = blood2016_tpm_ex_in$Row.names
blood2016_tpm_ex_in$Row.names = NULL

blood2016_tpm_ex_in$gene = rownames(blood2016_tpm_ex_in)


################
### cleanCd8 ###
################


cleanCd8_counts_gb = read.table("featureCounts_output/cleanCd8_counts_gb.txt", header = T, row.names = 1)
gb_length = cleanCd8_counts_gb$Length
cleanCd8_counts_gb = cleanCd8_counts_gb[,6:ncol(cleanCd8_counts_gb)]
coln_tmp = sub('.*CD', 'CD', colnames(cleanCd8_counts_gb))
coln_tmp = sub('_CKDL.*', '', coln_tmp)
targetFile = read.table('cleanCd8/targetFile.txt', header = T)
targetFile$cd = sub('.ReadsPerGene.out.tab.rawCounts', '', targetFile$files)
colnames(cleanCd8_counts_gb) = targetFile[match(coln_tmp, targetFile$cd),]$label
cleanCd8_counts_gb = cleanCd8_counts_gb[, c('TN.C.4wts.1', 'TN.C.4wts.2', 'TN.C.4wts.3', 
                                            'VM.C.4wts.1', 'VM.C.4wts.2', 'VM.C.4wts.3',
                                            'TN.D.4wts.1', 'TN.D.4wts.2', 'TN.D.4wts.3',
                                            'VM.D.4wts.1', 'VM.D.4wts.2', 'VM.D.4wts.3')]
colnames(cleanCd8_counts_gb) = c("Clean_TN_rep1", "Clean_TN_rep2", "Clean_TN_rep3",  
                                 "Clean_VM_rep1", "Clean_VM_rep2", "Clean_VM_rep3",
                                 "Dirty_TN_rep1", "Dirty_TN_rep2", "Dirty_TN_rep3", 
                                 "Dirty_VM_rep1", "Dirty_VM_rep2", "Dirty_VM_rep3")

cleanCd8_counts_exon = read.table("featureCounts_output/cleanCd8_counts_exon.txt", header = T, row.names = 1)
exon_length = cleanCd8_counts_exon$Length
cleanCd8_counts_exon = cleanCd8_counts_exon[,6:ncol(cleanCd8_counts_exon)]
coln_tmp = sub('.*CD', 'CD', colnames(cleanCd8_counts_exon))
coln_tmp = sub('_CKDL.*', '', coln_tmp)
targetFile = read.table('cleanCd8/targetFile.txt', header = T)
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


cleanCd8_counts_intron <- cleanCd8_counts_gb - cleanCd8_counts_exon
intron_length = gb_length - exon_length

cleanCd8_counts_intron_noNeg <- cleanCd8_counts_intron[rowSums(cleanCd8_counts_intron < 0) == 0,]
cleanCd8_counts_exon_noNeg <- cleanCd8_counts_exon[rowSums(cleanCd8_counts_intron < 0) == 0,]

intron_length_noNeg = intron_length[rowSums(cleanCd8_counts_intron < 0) == 0]
exon_length_noNeg = exon_length[rowSums(cleanCd8_counts_intron < 0) == 0]


cleanCd8_counts_intron_wIntron = cleanCd8_counts_intron_noNeg[intron_length_noNeg != 0, ]
cleanCd8_counts_exon_wIntron = cleanCd8_counts_exon_noNeg[intron_length_noNeg != 0, ]

intron_length_wIntron = intron_length_noNeg[intron_length_noNeg != 0]
exon_length_wIntron = exon_length_noNeg[intron_length_noNeg != 0]


cleanCd8_tpm_exon_wIntron <- apply(cleanCd8_counts_exon_wIntron, 2, function(x) tpm(x, exon_length_wIntron))
cleanCd8_tpm_intron_wIntron <- apply(cleanCd8_counts_intron_wIntron, 2, function(x) tpm(x, intron_length_wIntron))

cleanCd8_tpm_ex_in = merge(cleanCd8_tpm_exon_wIntron, cleanCd8_tpm_intron_wIntron, by = 0, suffixes = c(".exon",".intron"))
rownames(cleanCd8_tpm_ex_in) = cleanCd8_tpm_ex_in$Row.names
cleanCd8_tpm_ex_in$Row.names = NULL

cleanCd8_tpm_ex_in$gene = rownames(cleanCd8_tpm_ex_in)



################
### mir29cre ###
################

mir29cre_counts_gb = read.table("featureCounts_output/mir29cre_counts_gb.txt", header = T, row.names = 1)
gb_length = mir29cre_counts_gb$Length
mir29cre_counts_gb = mir29cre_counts_gb[,6:ncol(mir29cre_counts_gb)]
colnames(mir29cre_counts_gb) = sub("X.home.hz543.data.1087.Ciaran.mir29KO.RNAseq.redo_220406.hisat2.", "", colnames(mir29cre_counts_gb))
colnames(mir29cre_counts_gb) = c('1.WT.TN', '1.WT.MP', '1.mCre.TN', '1.mCre.MP', '1.pCre.TN', '1.pCre.MP',
                                 '2.WT.TN', '2.WT.MP', '2.mCre.TN', '2.mCre.MP', '2.pCre.TN', '2.pCre.MP')
mir29cre_counts_gb = mir29cre_counts_gb[, c('1.WT.TN', '1.WT.MP', '1.mCre.TN', '1.mCre.MP',
                                            '2.WT.TN', '2.WT.MP', '2.mCre.TN', '2.mCre.MP')]

mir29cre_counts_exon = read.table("featureCounts_output/mir29cre_counts_exon.txt", header = T, row.names = 1)
exon_length = mir29cre_counts_exon$Length
mir29cre_counts_exon = mir29cre_counts_exon[,6:ncol(mir29cre_counts_exon)]
colnames(mir29cre_counts_exon) = sub("X.home.hz543.data.1087.Ciaran.mir29KO.RNAseq.redo_220406.hisat2.", "", colnames(mir29cre_counts_exon))
colnames(mir29cre_counts_exon) = c('1.WT.TN', '1.WT.MP', '1.mCre.TN', '1.mCre.MP', '1.pCre.TN', '1.pCre.MP',
                                   '2.WT.TN', '2.WT.MP', '2.mCre.TN', '2.mCre.MP', '2.pCre.TN', '2.pCre.MP')
mir29cre_counts_exon = mir29cre_counts_exon[, c('1.WT.TN', '1.WT.MP', '1.mCre.TN', '1.mCre.MP',
                                                '2.WT.TN', '2.WT.MP', '2.mCre.TN', '2.mCre.MP')]


#table(rownames(mir29cre_counts_gb) == rownames(mir29cre_counts_exon))

mir29cre_counts_intron <- mir29cre_counts_gb - mir29cre_counts_exon
intron_length = gb_length - exon_length
# <0 intron counts: mir29cre_counts_intron[rowSums(mir29cre_counts_intron < 0) > 0,]
#write.table(mir29cre_counts_intron[rowSums(mir29cre_counts_intron < 0) > 0,], 'test.txt', quote = F)'

mir29cre_counts_intron_noNeg <- mir29cre_counts_intron[rowSums(mir29cre_counts_intron < 0) == 0,]
mir29cre_counts_exon_noNeg <- mir29cre_counts_exon[rowSums(mir29cre_counts_intron < 0) == 0,]

intron_length_noNeg = intron_length[rowSums(mir29cre_counts_intron < 0) == 0]
exon_length_noNeg = exon_length[rowSums(mir29cre_counts_intron < 0) == 0]


mir29cre_counts_intron_wIntron = mir29cre_counts_intron_noNeg[intron_length_noNeg != 0, ]
mir29cre_counts_exon_wIntron = mir29cre_counts_exon_noNeg[intron_length_noNeg != 0, ]

intron_length_wIntron = intron_length_noNeg[intron_length_noNeg != 0]
exon_length_wIntron = exon_length_noNeg[intron_length_noNeg != 0]


mir29cre_tpm_exon_wIntron <- apply(mir29cre_counts_exon_wIntron, 2, function(x) tpm(x, exon_length_wIntron))
mir29cre_tpm_intron_wIntron <- apply(mir29cre_counts_intron_wIntron, 2, function(x) tpm(x, intron_length_wIntron))

mir29cre_tpm_ex_in = merge(mir29cre_tpm_exon_wIntron, mir29cre_tpm_intron_wIntron, by = 0, suffixes = c(".exon",".intron"))
rownames(mir29cre_tpm_ex_in) = mir29cre_tpm_ex_in$Row.names
mir29cre_tpm_ex_in$Row.names = NULL

mir29cre_tpm_ex_in$gene = rownames(mir29cre_tpm_ex_in)


##############
### cd8dev ###
##############

library("readxl")
filenames <- read_excel('cd8dev/Filenames.xlsx', col_names = F)
filenames$DGnum = gsub('N.ReadsPerGene.out.tab.rawCounts', '', filenames$...3)
filenames$DGnum = gsub('c', '', filenames$DGnum)

cd8dev_counts_gb = read.table("featureCounts_output/cd8dev_counts_gb.txt", header = T, row.names = 1)
gb_length = cd8dev_counts_gb$Length
cd8dev_counts_gb = cd8dev_counts_gb[,6:ncol(cd8dev_counts_gb)]
colnames(cd8dev_counts_gb) = gsub('.*.DG', 'DG', colnames(cd8dev_counts_gb))
colnames(cd8dev_counts_gb) = gsub('N.*', '', colnames(cd8dev_counts_gb))
#table(as.data.frame(filenames[match(colnames(cd8dev_counts_gb),filenames$DGnum),'DGnum'])$DGnum == colnames(cd8dev_counts_gb))
colnames(cd8dev_counts_gb) = as.data.frame(filenames[match(colnames(cd8dev_counts_gb),filenames$DGnum),'...2'])$'...2'
cd8dev_counts_gb = cd8dev_counts_gb[, c("ACD81", "NCD81", 
                                        "ACD82", "NCD82", 
                                        "ACD83", "NCD83", 
                                        "ACD84", "NCD84")]

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


#table(rownames(cd8dev_counts_gb) == rownames(cd8dev_counts_exon))

cd8dev_counts_intron <- cd8dev_counts_gb - cd8dev_counts_exon
intron_length = gb_length - exon_length

cd8dev_counts_intron_noNeg <- cd8dev_counts_intron[rowSums(cd8dev_counts_intron < 0) == 0,]
cd8dev_counts_exon_noNeg <- cd8dev_counts_exon[rowSums(cd8dev_counts_intron < 0) == 0,]

intron_length_noNeg = intron_length[rowSums(cd8dev_counts_intron < 0) == 0]
exon_length_noNeg = exon_length[rowSums(cd8dev_counts_intron < 0) == 0]


cd8dev_counts_intron_wIntron = cd8dev_counts_intron_noNeg[intron_length_noNeg != 0, ]
cd8dev_counts_exon_wIntron = cd8dev_counts_exon_noNeg[intron_length_noNeg != 0, ]

intron_length_wIntron = intron_length_noNeg[intron_length_noNeg != 0]
exon_length_wIntron = exon_length_noNeg[intron_length_noNeg != 0]


cd8dev_tpm_exon_wIntron <- apply(cd8dev_counts_exon_wIntron, 2, function(x) tpm(x, exon_length_wIntron))
cd8dev_tpm_intron_wIntron <- apply(cd8dev_counts_intron_wIntron, 2, function(x) tpm(x, intron_length_wIntron))

cd8dev_tpm_ex_in = merge(cd8dev_tpm_exon_wIntron, cd8dev_tpm_intron_wIntron, by = 0, suffixes = c(".exon",".intron"))
rownames(cd8dev_tpm_ex_in) = cd8dev_tpm_ex_in$Row.names
cd8dev_tpm_ex_in$Row.names = NULL

cd8dev_tpm_ex_in$gene = rownames(cd8dev_tpm_ex_in)


################## datasets loaded ##################



# combining data sets

tpm_no0_all = Reduce(function(x,y) merge(x = x, y = y, by = 'gene'), 
                     list(mir29_tpm_ex_in, erin2015_tpm_ex_in, blood2016_tpm_ex_in, cleanCd8_tpm_ex_in, mir29cre_tpm_ex_in, cd8dev_tpm_ex_in)) 

row.names(tpm_no0_all) <- tpm_no0_all$gene
tpm_no0_all[,'gene'] <- NULL

tpm_no0_all_exon = tpm_no0_all[ , grepl( ".exon" , names( tpm_no0_all ) ) ]
tpm_no0_all_intron = tpm_no0_all[ , grepl( ".intron" , names( tpm_no0_all ) ) ]


# initialize decay rate matrix 
decay_rates = data.frame(matrix(nrow = nrow(tpm_no0_all_exon), ncol = ncol(tpm_no0_all_exon)))
colnames(decay_rates) = gsub('.exon', '', colnames(tpm_no0_all_exon))
rownames(decay_rates) = rownames(tpm_no0_all_exon)

for (i in 1:nrow(decay_rates)){
  for (j in 1:ncol(decay_rates)){
    decay_rates[i,j] = tpm_no0_all_intron[i,j] / tpm_no0_all_exon[i,j]
  }
}
write.table(tpm_no0_all, paste0(outputdir, "tpm_no0_all.txt"), quote = F, sep = "\t")
write.table(decay_rates, paste0(outputdir, "decay_rates.txt"), quote = F, sep = "\t")



# matching to smRNAseq, to be used as input for miR-Inf 
mirna_exp = read.table('inputs/mirna_exp_matrix/counts_vst_blindT_combat_expressed1000.txt') # this was generated in Fig3 code

decayRates_mtch = data.frame(matrix(nrow = nrow(decay_rates), ncol = 41)) #58 64
#colnames(decayRates_mtch)[1:2] = c('ilin28b_spleen', 'wt_spleen')
rownames(decayRates_mtch) = rownames(decay_rates)
#decayRates_mtch$ilin28b_spleen = rowMeans(decay_rates[,c("Lin28b_TNa", "Lin28b_MPa", "Lin28b_TNb", "Lin28b_MPb", "Lin28b_TNc", "Lin28b_MPc", "Lin28b_TNd", "Lin28b_MPd")])
#decayRates_mtch$wt_spleen = rowMeans(decay_rates[,c("WT_TNa", "WT_MPa", "WT_TNb", "WT_MPb", "WT_TNc", "WT_MPc", "WT_TNd", "WT_MPd")])


# get mean across replicates to match smRNAseq data with different replicates
colnames(decayRates_mtch)[1:7] = c('adult_vm_1', 'adult_vm_2', 'adult_tn_1', 'adult_tn_2', 'neo_vm_1', 'neo_vm_2', 'neo_tn_1')
decayRates_mtch$adult_vm_1 = rowMeans(decay_rates[,c('Adult-1-VM', 'Adult-2-VM', 'Adult-3-VM')])
decayRates_mtch$adult_vm_2 = rowMeans(decay_rates[,c('Adult-1-VM', 'Adult-2-VM', 'Adult-3-VM')])
decayRates_mtch$adult_tn_1 = rowMeans(decay_rates[,c('Adult-1-TN', 'Adult-2-TN', 'Adult-3-TN')])
decayRates_mtch$adult_tn_2 = rowMeans(decay_rates[,c('Adult-1-TN', 'Adult-2-TN', 'Adult-3-TN')])
decayRates_mtch$neo_vm_1 = rowMeans(decay_rates[,c('NeoNC-1-VM', 'NeoNC-2-VM', 'NeoNC-3-VM')])
decayRates_mtch$neo_vm_2 = rowMeans(decay_rates[,c('NeoNC-1-VM', 'NeoNC-2-VM', 'NeoNC-3-VM')])
decayRates_mtch$neo_tn_1 = rowMeans(decay_rates[,c('NeoNC-1-TN', 'NeoNC-2-TN', 'NeoNC-3-TN')])

colnames(decayRates_mtch)[8:11] = c('adult_naive_ilm', 'adult_naive_neb', 'neo_naive_ilm', 'neo_naive_neb')
decayRates_mtch[,8:11] =  decay_rates[,c('adult_naive', 'adult_naive', 'neonate_naive', 'neonate_naive')]

colnames(decayRates_mtch)[12:13] = c('splenic_adult', 'splenic_neonate') 
decayRates_mtch[,12:13] =  decay_rates[,c('adult_spleen', 'neonatal_spleen')] 

colnames(decayRates_mtch)[14:25] = c("Clean_TN_rep1", "Clean_TN_rep2", "Clean_TN_rep3",  "Clean_VM_rep1", "Clean_VM_rep2", "Clean_VM_rep3","Dirty_TN_rep1", "Dirty_TN_rep2", "Dirty_TN_rep3", "Dirty_VM_rep1", "Dirty_VM_rep2", "Dirty_VM_rep3")
decayRates_mtch[,14:25] =  decay_rates[,c("Clean_TN_rep1", "Clean_TN_rep2", "Clean_TN_rep3",  "Clean_VM_rep1", "Clean_VM_rep2", "Clean_VM_rep3","Dirty_TN_rep1", "Dirty_TN_rep2", "Dirty_TN_rep3", "Dirty_VM_rep1", "Dirty_VM_rep2", "Dirty_VM_rep3")]

colnames(decayRates_mtch)[26:33] = c('mm1', 'mm2', 'mt1', 'mt2',
                                     'wm1', 'wm2', 'wt1', 'wt2')
decayRates_mtch[,26:33] =  decay_rates[,c('1.mCre.MP', '2.mCre.MP', '1.mCre.TN', '2.mCre.TN', 
                                          '1.WT.MP'  , '2.WT.MP'  , '1.WT.TN'  , '2.WT.TN'  )]

colnames(decayRates_mtch)[34:41] = c("ACD81", "NCD81", "ACD82", "NCD82", "ACD83", "NCD83", "ACD84", "NCD84")
decayRates_mtch[,34:41] =  decay_rates[,c("ACD81", "NCD81", "ACD82", "NCD82", "ACD83", "NCD83", "ACD84", "NCD84")]



decayRates_mtch_noNA <- na.omit(decayRates_mtch)
decayRates_mtch_noNA_noInf <- decayRates_mtch_noNA[!is.infinite(rowSums(decayRates_mtch_noNA)),]
decayRates_mtch_noNA_noInf_noall0 = decayRates_mtch_noNA_noInf[rowSums(decayRates_mtch_noNA_noInf[])>0,]


### additional filtering on lowly expressed genes
outputdir = 'inputs/exonIntron_matrix/'

tpm_all = read.table('inputs/exonIntron_matrix/tpm_no0_all.txt')
tpm_exon = tpm_all[ , grepl( ".exon" , names( tpm_all ) ) ]
tpm_exon_10exon = tpm_exon[rowSums(tpm_exon > 10) >= 1, ] # 18397 genes -> 7173 genes
shared_tmp = intersect(rownames(tpm_exon_10exon), rownames(decayRates_mtch_noNA_noInf_noall0))
decayRates_mtch_noNA_noInf_noall0_10exon = decayRates_mtch_noNA_noInf_noall0[shared_tmp, ]
write.table(decayRates_mtch_noNA_noInf_noall0_10exon, paste0(outputdir, "decay_rates_match_noNA_noInf_noall0_10exon.txt"), quote = F, sep = "\t")

targGene_decay_rates_match_noNA_noInf_noall0_10exon = rownames(decayRates_mtch_noNA_noInf_noall0_10exon) # 6323 genes
write.table(targGene_decay_rates_match_noNA_noInf_noall0_10exon, paste0(outputdir, "targGene_match_noNA_noInf_noall0_10exon.txt"), row.names = FALSE, col.names=FALSE, sep = "\t", quote = F)

