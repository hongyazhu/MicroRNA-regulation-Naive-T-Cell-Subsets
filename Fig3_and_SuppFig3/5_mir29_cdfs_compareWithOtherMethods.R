library(VennDiagram)
library(edgeR)
library(reshape2)
library(ggplot2)
library(tools)
library(readxl)

# melt message suppressed

####################
### read data in ###
####################

# from Targetscan, find miRNAs' seed sequence (used to find targets) 
mir_family_info <- read.table("targetscan/mouse/miR_Family_Info.txt", head = T, fill = T, sep = "\t") 
hsa_mir_family_info <- mir_family_info[substr(mir_family_info[,'MiRBase.ID'] ,1 ,3) == "mmu", ] 
hsa_con_mir_family_info <- hsa_mir_family_info[hsa_mir_family_info[,'Family.Conservation.'] >= 2, ] 
seed_con <- hsa_con_mir_family_info[,c('Seed.m8', 'MiRBase.ID')] 
mirna_seeds_con <- aggregate(seed_con[,2], by = list(seed = seed_con$Seed), FUN = paste) 
colnames(mirna_seeds_con) <- c("seed", "miRNA") 
mirna_seeds_con$miRNA_name = c('miR-106ab,17,20ab,93-5p,6383','miR-141,200a-3p','miR-132,212-3p','miR-451a','miR-490-3p','miR-191-5p','miR-124-3p.1','miR-18ab-5p','miR-291a,294,295,302abd-3p','miR-200bc,429-3p','miR-216a-5p','miR-216b-5p','miR-365-3p','miR-101a-3p.1','miR-144-3p','miR-181abcd-5p','miR-140-3p.2,miR-497b','miR-100,99ab-5p','miR-10ab-5p','miR-217-5p','miR-193ab-3p','miR-383-5p.2','miR-29abc-3p','miR-15ab,16,1907,195a,322,497a-5p,miR-195b,6342,6353,6419','miR-129-1,2-3p','miR-22-3p','miR-21a,590-5p,miR-21c','miR-196ab-5p','miR-130ab,301ab-3p,miR-130c,6341,6389,721','miR-302c-3p','miR-140-5p,miR-876-3p','miR-142a-5p','miR-425-5p,miR-489-3p','miR-183-5p','miR-135ab-5p','miR-455-5p,miR-5129-3p','miR-25,363,367,92ab-3p,miR-32-5p','miR-128-3p,miR-6539','miR-802-5p','miR-199ab-3p','miR-455-3p.1','miR-148ab,152-3p','miR-140-3p.1','miR-338-3p','miR-199ab-5p','miR-125ab,351-5p,miR-6367,6394','miR-205-5p','miR-212-5p','miR-551b-3p','miR-126a-3p.1','miR-187-3p','miR-139-5p','miR-150-5p,miR-5127','miR-9-5p','miR-203-3p.1','miR-146ab-5p','miR-143-3p','let-7abcdefgik,98-5p,miR-1961','miR-190ab-5p','miR-383-5p.1','miR-219a-5p','miR-103,107-3p','miR-214-5p','miR-221,222-3p,miR-1928','miR-138-5p','miR-7ab-5p','miR-1a,206-3p,miR-1957b,6349,6382','miR-184-3p','miR-122-5p','miR-31-5p','miR-34abc,449ac-5p,-miR-449b','miR-24-3p,miR-5124b,6361,6369,6410,6413','miR-193a-5p','miR-30abcde,384-5p','miR-194-5p','miR-126a-3p.2','miR-142a-3p.1','miR-223-3p','miR-19ab-3p','miR-208ab-3p','miR-499-5p','miR-124,5624-3p,miR-6540-5p','miR-155-5p','miR-101ab-3p','miR-142a-3p.2','miR-137-3p','miR-26ab-5p','miR-27ab-3p','miR-23ab-3p','miR-145a-5p,miR-145b','miR-204,211-5p,miR-7670-3p','miR-202-5p','miR-203-3p.2','miR-192,215-5p','miR-455-3p.2,miR-682','miR-153-3p','miR-33-5p','miR-183-5p.2','miR-133a-3p.1','miR-147-3p','miR-210-3p','miR-218-5p,miR-7002-3p','miR-182-5p','miR-96-5p','miR-133ab-3p,miR-133c','miR-375-3p','miR-129-5p')

# Targetscan
ts_nopct = read.table('inputs/prior/ts_ct02_b.tsv')
ts = read.table('inputs/prior/ts_ct02_pct75_b.tsv')

### find weak targets (to be excluded in the analysis) 
weak_targets <- read.table("/home/hz543/data/targetscan/mouse/Targetscan_mouse_weak_targets_genes_names.txt") 
colnames(weak_targets) <- c("transcript_id_ver", "gene_name","seed", "mirna_name") 
weak_targets <- weak_targets[, c("gene_name","seed")] 
weak_targets <- weak_targets[weak_targets$seed %in% colnames(ts_nopct), ]

# miRDB - downloaded from miRDB website
miRDB_mir29a = read_excel('miRDB_mir29_predictions.xlsx', sheet = 'mmu-miR-29a-3p')
miRDB_mir29b = read_excel('miRDB_mir29_predictions.xlsx', sheet = 'mmu-miR-29b-3p')
miRDB_mir29c = read_excel('miRDB_mir29_predictions.xlsx', sheet = 'mmu-miR-29c-3p')

miRDB_mir29abc = unique(c(miRDB_mir29a$`Gene Symbol`, miRDB_mir29b$`Gene Symbol`, miRDB_mir29c$`Gene Symbol`))

# miRWalk - downloaded from miRWalk website
miRWalk_mir29a = read.csv('miRWalk_miR-29a-3p_Targets.csv', sep = ',', header = T)
miRWalk_mir29b = read.csv('miRWalk_miR-29b-3p_Targets.csv', sep = ',', header = T)
miRWalk_mir29c = read.csv('miRWalk_miR-29c-3p_Targets.csv', sep = ',', header = T)

miRWalk_mir29abc = unique(c(miRWalk_mir29a$genesymbol, miRWalk_mir29b$genesymbol, miRWalk_mir29c$genesymbol))

# anamiR - generated in run_anamiR
anamiR_res = read.table('anamiR/res.txt', sep = '\t', header = T)
anamiR_mir29 = unique(anamiR_res[anamiR_res$miRNA_name == 'miR-29abc-3p','Gene'])

#########################
### get relevant data ### 
#########################

library("readxl")
filenames <- read_excel('Filenames.xlsx', col_names = F)
filenames$DGnum = gsub('N.ReadsPerGene.out.tab.rawCounts', '', filenames$...3)
filenames$DGnum = gsub('c', '', filenames$DGnum)
cd8dev_counts_wholegene = read.table("cd8dev/counts.txt", header = T, row.names = 1)
wholegene_length = cd8dev_counts_wholegene$Length
cd8dev_counts_wholegene = cd8dev_counts_wholegene[,6:ncol(cd8dev_counts_wholegene)]
colnames(cd8dev_counts_wholegene) = gsub('.*.DG', 'DG', colnames(cd8dev_counts_wholegene))
colnames(cd8dev_counts_wholegene) = gsub('N.*', '', colnames(cd8dev_counts_wholegene))
#table(as.data.frame(filenames[match(colnames(cd8dev_counts_wholegene),filenames$DGnum),'DGnum'])$DGnum == colnames(cd8dev_counts_wholegene))
colnames(cd8dev_counts_wholegene) = as.data.frame(filenames[match(colnames(cd8dev_counts_wholegene),filenames$DGnum),'...2'])$'...2'
cd8dev_counts_wholegene_rel = cd8dev_counts_wholegene[,c("ACD81", "NCD81", "NCD82", "NCD83", "NCD84", "ACD82",  "ACD83", "ACD84")]

cd8dev_counts_wholegene_rel_cpm = data.frame(cpm(cd8dev_counts_wholegene_rel))
cd8dev_counts_wholegene_rel_cpm_noLow = cd8dev_counts_wholegene_rel_cpm[rowSums(cd8dev_counts_wholegene_rel_cpm < 1) == 0 , , drop = FALSE] # 12874 genes
cd8dev_counts_wholegene_rel_cpm_noLow = rownames(cd8dev_counts_wholegene_rel_cpm_noLow)

plot_barplot_lineplots_cd8devproject <- function(matrix, seed, network_path, plot_title, file_name){
  
  decay_rates = matrix # it's called decay_rates in this code, but it can be decay_rates, wholegene_cpm to be plotted
  
  ######################
  ### set parameters ### 
  ######################
  
  
  miRNA_name = mirna_seeds_con[mirna_seeds_con$seed == seed, 'miRNA_name']
  
  
  dir.create(file.path(paste0('make_figures/ct02_ms2_4ormore/benchmark/mir29/', 
                              "/cd8devproject_validatetargets_pct/", miRNA_name, "/")), recursive = TRUE, showWarnings = FALSE)
  
  ###################
  ### get targets ### 
  ###################
  
  # weak targets from targetscan - to be excluded in the analysis 
  if (!grepl('\\|', seed)){
    miRNA_weak_ts <- weak_targets[weak_targets$seed == seed,]$gene_name
  } else {
    seed1 = strsplit(seed,split='|', fixed=T)[[1]][1]
    seed2 = strsplit(seed,split='|', fixed=T)[[1]][2]
    miRNA_weak_ts1 <- weak_targets[weak_targets$seed == seed1,]$gene_name
    miRNA_weak_ts2 <- weak_targets[weak_targets$seed == seed2,]$gene_name
    miRNA_weak_ts = union(miRNA_weak_ts1, miRNA_weak_ts2)
  }
  #miRNA_weak_ts <- weak_targets_geneID[weak_targets_geneID$seed == seed,]$gene_id
  
  
  
  ### Targetscan targets
  # context++ score < –0.2, pct < 0.75
  if (!grepl('\\|', seed)){
    miRNA_ts <- rownames(ts)[ts[, seed] == 1]
  } else {
    seed1 = strsplit(seed,split='|', fixed=T)[[1]][1]
    seed2 = strsplit(seed,split='|', fixed=T)[[1]][2]
    miRNA_ts1 <- rownames(ts)[ts[, seed1] == 1]
    miRNA_ts2 <- rownames(ts)[ts[, seed2] == 1]
    miRNA_ts = union(miRNA_ts1, miRNA_ts2)
  }
  # context++ score < –0.2
  if (!grepl('\\|', seed)){
    miRNA_tsnopct <- rownames(ts_nopct)[ts_nopct[, seed] == 1]
  } else {
    seed1 = strsplit(seed,split='|', fixed=T)[[1]][1]
    seed2 = strsplit(seed,split='|', fixed=T)[[1]][2]
    miRNA_ts1nopct <- rownames(ts_nopct)[ts_nopct[, seed1] == 1]
    miRNA_ts2nopct <- rownames(ts_nopct)[ts_nopct[, seed2] == 1]
    miRNA_tsnopct = union(miRNA_ts1nopct, miRNA_ts2nopct)
  }
  
  
  # include only genes that are not lowly expressed
  miRNA_ts_notLow = intersect(miRNA_ts, cd8dev_counts_wholegene_rel_cpm_noLow)
  miRNA_tsnopct_notLow = intersect(miRNA_tsnopct, cd8dev_counts_wholegene_rel_cpm_noLow)
  
  # remove weak targets
  # miRNA_ts_noWeak = setdiff(miRNA_ts_notLow, miRNA_weak_ts)
  
  stopifnot(miRNA_name == 'miR-29abc-3p')
  
  # miRDB
  miRNA_miRDB_mir29abc_notLow = intersect(miRDB_mir29abc, cd8dev_counts_wholegene_rel_cpm_noLow)
  miRNA_miRDB_mir29abc_alsoTs = intersect(miRNA_miRDB_mir29abc_notLow, miRNA_ts_notLow)
  
  # miRWalk
  miRNA_miRWalk_mir29abc_notLow = intersect(miRWalk_mir29abc, cd8dev_counts_wholegene_rel_cpm_noLow)
  miRNA_miRWalk_mir29abc_alsoTs = intersect(miRNA_miRWalk_mir29abc_notLow, miRNA_ts_notLow)
  
  # anamiR
  miRNA_anamiR_mir29abc_notLow = intersect(anamiR_mir29, cd8dev_counts_wholegene_rel_cpm_noLow)
  miRNA_anamiR_mir29abc_alsoTs = intersect(miRNA_anamiR_mir29abc_notLow, miRNA_ts_notLow)
  
  
  ### Inferelator targets
  network = read.table(network_path, sep = "\t", header = T)
  miRNA_inf = network[network$TF == seed, 'Target'] 
  
  # include only genes that are not lowly expressed
  miRNA_inf_notLow = intersect(miRNA_inf, cd8dev_counts_wholegene_rel_cpm_noLow)
  
  # remove weak targets
  miRNA_inf_noWeak = setdiff(miRNA_inf_notLow, miRNA_weak_ts)
  
  
  ### intersect of Inferelator and Targetscan targets
  miRNA_both = intersect(miRNA_inf_noWeak, miRNA_ts_notLow)
  miRNA_bothnopct = intersect(miRNA_inf_noWeak, miRNA_tsnopct_notLow)
  
  
  ###################################
  ### get decay rates for targets ### 
  ###################################
  
  # decay rates relevant to this miRNA
  decay_rates_rel = decay_rates[,c("ACD81", "NCD81", "NCD82", "NCD83", "NCD84", "ACD82",  "ACD83", "ACD84")]
  
  # remove genes with NAs and Inf decay rates (dr)
  decay_rates_rel_noNA <- na.omit(decay_rates_rel)
  decay_rates_rel_noNA_noInf <- decay_rates_rel_noNA[!is.infinite(rowSums(decay_rates_rel_noNA)),]
  
  ### select decay rates with targets
  miRNA_ts_dr = decay_rates_rel_noNA_noInf[intersect(miRNA_ts_notLow, rownames(decay_rates_rel_noNA_noInf)), ]
  miRNA_tsnopct_dr = decay_rates_rel_noNA_noInf[intersect(miRNA_tsnopct_notLow, rownames(decay_rates_rel_noNA_noInf)), ]
  miRNA_inf_dr = decay_rates_rel_noNA_noInf[intersect(miRNA_inf_noWeak, rownames(decay_rates_rel_noNA_noInf)), ] # not plotting as it completely overlap with both
  miRNA_both_dr = decay_rates_rel_noNA_noInf[intersect(miRNA_both, rownames(decay_rates_rel_noNA_noInf)), ]
  miRNA_bothnopct_dr = decay_rates_rel_noNA_noInf[intersect(miRNA_bothnopct, rownames(decay_rates_rel_noNA_noInf)), ]
  miRNA_miRDB_dr = decay_rates_rel_noNA_noInf[intersect(miRNA_miRDB_mir29abc_notLow, rownames(decay_rates_rel_noNA_noInf)), ]
  miRNA_miRDBTs_dr = decay_rates_rel_noNA_noInf[intersect(miRNA_miRDB_mir29abc_alsoTs, rownames(decay_rates_rel_noNA_noInf)), ]
  miRNA_miRWalk_dr = decay_rates_rel_noNA_noInf[intersect(miRNA_miRWalk_mir29abc_notLow, rownames(decay_rates_rel_noNA_noInf)), ]
  miRNA_miRWalkTs_dr = decay_rates_rel_noNA_noInf[intersect(miRNA_miRWalk_mir29abc_alsoTs, rownames(decay_rates_rel_noNA_noInf)), ]
  miRNA_anamiR_dr = decay_rates_rel_noNA_noInf[intersect(miRNA_anamiR_mir29abc_notLow, rownames(decay_rates_rel_noNA_noInf)), ]
  miRNA_anamiRTs_dr = decay_rates_rel_noNA_noInf[intersect(miRNA_anamiR_mir29abc_alsoTs, rownames(decay_rates_rel_noNA_noInf)), ]
  
  
  ### select decay rates with background genes
  # include only genes that are not lowly expressed
  bg_notLow = intersect(rownames(decay_rates_rel_noNA_noInf), cd8dev_counts_wholegene_rel_cpm_noLow)
  
  # remove weak targets
  bg_noWeak = setdiff(bg_notLow, miRNA_weak_ts)
  
  # remove all targets
  bg_noTargets = setdiff(bg_noWeak, miRNA_ts_notLow)
  bg_noTargets = setdiff(bg_noTargets, miRNA_tsnopct_notLow)
  bg_noTargets = setdiff(bg_noTargets, miRNA_both)
  bg_noTargets = setdiff(bg_noTargets, miRNA_bothnopct)
  bg_noTargets = setdiff(bg_noTargets, miRNA_miRDB_mir29abc_notLow)
  bg_noTargets = setdiff(bg_noTargets, miRNA_miRWalk_mir29abc_notLow) # removes a lot of genes
  bg_noTargets = setdiff(bg_noTargets, miRNA_anamiR_mir29abc_notLow)
  
  miRNA_bg_dr = decay_rates_rel_noNA_noInf[intersect(bg_noTargets, rownames(decay_rates_rel_noNA_noInf)), ]
  
  miRNA_ts_dr_toMelt = data.frame(log10(miRNA_ts_dr))
  miRNA_tsnopct_dr_toMelt = data.frame(log10(miRNA_tsnopct_dr))
  miRNA_inf_dr_toMelt = data.frame(log10(miRNA_inf_dr))
  miRNA_both_dr_toMelt = data.frame(log10(miRNA_both_dr))
  miRNA_bothnopct_dr_toMelt = data.frame(log10(miRNA_bothnopct_dr))
  miRNA_miRDB_dr_toMelt = data.frame(log10(miRNA_miRDB_dr))
  miRNA_miRDBTs_dr_toMelt = data.frame(log10(miRNA_miRDBTs_dr))
  miRNA_miRWalk_dr_toMelt = data.frame(log10(miRNA_miRWalk_dr))
  miRNA_miRWalkTs_dr_toMelt = data.frame(log10(miRNA_miRWalkTs_dr))
  miRNA_anamiR_dr_toMelt = data.frame(log10(miRNA_anamiR_dr))
  miRNA_anamiRTs_dr_toMelt = data.frame(log10(miRNA_anamiRTs_dr))
  miRNA_bg_dr_toMelt = data.frame(log10(miRNA_bg_dr))
  
  miRNA_ts_dr_toMelt$method = 'Targetscan(conserved)'
  miRNA_tsnopct_dr_toMelt$method = 'Targetscan'
  miRNA_inf_dr_toMelt$method = 'Inferelator'
  miRNA_both_dr_toMelt$method = 'Both(conserved)'
  miRNA_bothnopct_dr_toMelt$method = 'Both'
  miRNA_miRDB_dr_toMelt$method = 'miRDB'
  miRNA_miRDBTs_dr_toMelt$method = 'miRBD+Targetscan'
  miRNA_miRWalk_dr_toMelt$method = 'miRWalk'
  miRNA_miRWalkTs_dr_toMelt$method = 'miRWalk+Targetscan'
  miRNA_anamiR_dr_toMelt$method = 'anamiR'
  miRNA_anamiRTs_dr_toMelt$method = 'anamiR+Targetscan'
  miRNA_bg_dr_toMelt$method = 'Background'
  
  suppressMessages({
    miRNA_ts_dr_melt = melt(miRNA_ts_dr_toMelt)
    miRNA_tsnopct_dr_melt = melt(miRNA_tsnopct_dr_toMelt)
    miRNA_inf_dr_melt = melt(miRNA_inf_dr_toMelt)
    miRNA_both_dr_melt = melt(miRNA_both_dr_toMelt)
    miRNA_bothnopct_dr_melt = melt(miRNA_bothnopct_dr_toMelt)
    miRNA_miRDB_dr_melt = melt(miRNA_miRDB_dr_toMelt)
    miRNA_miRDBTs_dr_melt = melt(miRNA_miRDBTs_dr_toMelt)
    miRNA_miRWalk_dr_melt = melt(miRNA_miRWalk_dr_toMelt)
    miRNA_miRWalkTs_dr_melt = melt(miRNA_miRWalkTs_dr_toMelt)
    miRNA_anamiR_dr_melt = melt(miRNA_anamiR_dr_toMelt)
    miRNA_anamiRTs_dr_melt = melt(miRNA_anamiRTs_dr_toMelt)
    miRNA_bg_dr_melt = melt(miRNA_bg_dr_toMelt)
  })
  
  miRNA_dr_melt = rbind(miRNA_ts_dr_melt, miRNA_tsnopct_dr_melt, miRNA_inf_dr_melt, miRNA_both_dr_melt, miRNA_bothnopct_dr_melt, 
                        miRNA_miRDB_dr_melt, miRNA_miRDBTs_dr_melt, miRNA_miRWalk_dr_melt, miRNA_miRWalkTs_dr_melt, miRNA_anamiR_dr_melt, miRNA_anamiRTs_dr_melt,
                        miRNA_bg_dr_melt)
  miRNA_dr_melt_finite <- miRNA_dr_melt[!is.infinite(miRNA_dr_melt$value),]
  
  ts_ntargets = nrow(miRNA_ts_dr)
  tsnopct_ntargets = nrow(miRNA_tsnopct_dr)
  inf_ntargets = nrow(miRNA_inf_dr)
  both_ntargets = nrow(miRNA_both_dr)
  bothnopct_ntargets = nrow(miRNA_bothnopct_dr)
  miRDB_ntargets = nrow(miRNA_miRDB_dr)
  miRDBTs_ntargets = nrow(miRNA_miRDBTs_dr)
  miRWalk_ntargets = nrow(miRNA_miRWalk_dr)
  miRWalkTs_ntargets = nrow(miRNA_miRWalkTs_dr)
  anamiR_ntargets = nrow(miRNA_anamiR_dr)
  anamiRTs_ntargets = nrow(miRNA_anamiRTs_dr)
  bg_ntargets = nrow(miRNA_bg_dr)
  
  
  # plot fold changes
  miRNA_inf_dr_fc = miRNA_inf_dr
  miRNA_ts_dr_fc = miRNA_ts_dr
  miRNA_tsnopct_dr_fc = miRNA_tsnopct_dr
  miRNA_both_dr_fc = miRNA_both_dr
  miRNA_bothnopct_dr_fc = miRNA_bothnopct_dr
  miRNA_miRDB_dr_fc = miRNA_miRDB_dr
  miRNA_miRDBTs_dr_fc = miRNA_miRDBTs_dr
  miRNA_miRWalk_dr_fc = miRNA_miRWalk_dr
  miRNA_miRWalkTs_dr_fc = miRNA_miRWalkTs_dr
  miRNA_anamiR_dr_fc = miRNA_anamiR_dr
  miRNA_anamiRTs_dr_fc = miRNA_anamiRTs_dr
  miRNA_bg_dr_fc = miRNA_bg_dr
  
  miRNA_inf_dr_fc$fc = log2(rowMeans(miRNA_inf_dr_fc[,c("NCD81", "NCD82", "NCD83")])/rowMeans(miRNA_inf_dr_fc[,c("ACD81", "ACD82", "ACD83")]))
  miRNA_ts_dr_fc$fc = log2(rowMeans(miRNA_ts_dr_fc[,c("NCD81", "NCD82", "NCD83")])/rowMeans(miRNA_ts_dr_fc[,c("ACD81", "ACD82", "ACD83")]))
  miRNA_tsnopct_dr_fc$fc = log2(rowMeans(miRNA_tsnopct_dr_fc[,c("NCD81", "NCD82", "NCD83")])/rowMeans(miRNA_tsnopct_dr_fc[,c("ACD81", "ACD82", "ACD83")]))
  miRNA_both_dr_fc$fc = log2(rowMeans(miRNA_both_dr_fc[,c("NCD81", "NCD82", "NCD83")])/rowMeans(miRNA_both_dr_fc[,c("ACD81", "ACD82", "ACD83")]))
  miRNA_bothnopct_dr_fc$fc = log2(rowMeans(miRNA_bothnopct_dr_fc[,c("NCD81", "NCD82", "NCD83")])/rowMeans(miRNA_bothnopct_dr_fc[,c("ACD81", "ACD82", "ACD83")]))
  miRNA_miRDB_dr_fc$fc = log2(rowMeans(miRNA_miRDB_dr_fc[,c("NCD81", "NCD82", "NCD83")])/rowMeans(miRNA_miRDB_dr_fc[,c("ACD81", "ACD82", "ACD83")]))
  miRNA_miRDBTs_dr_fc$fc = log2(rowMeans(miRNA_miRDBTs_dr_fc[,c("NCD81", "NCD82", "NCD83")])/rowMeans(miRNA_miRDBTs_dr_fc[,c("ACD81", "ACD82", "ACD83")]))
  miRNA_miRWalk_dr_fc$fc = log2(rowMeans(miRNA_miRWalk_dr_fc[,c("NCD81", "NCD82", "NCD83")])/rowMeans(miRNA_miRWalk_dr_fc[,c("ACD81", "ACD82", "ACD83")]))
  miRNA_miRWalkTs_dr_fc$fc = log2(rowMeans(miRNA_miRWalkTs_dr_fc[,c("NCD81", "NCD82", "NCD83")])/rowMeans(miRNA_miRWalkTs_dr_fc[,c("ACD81", "ACD82", "ACD83")]))
  miRNA_anamiR_dr_fc$fc = log2(rowMeans(miRNA_anamiR_dr_fc[,c("NCD81", "NCD82", "NCD83")])/rowMeans(miRNA_anamiR_dr_fc[,c("ACD81", "ACD82", "ACD83")]))
  miRNA_anamiRTs_dr_fc$fc = log2(rowMeans(miRNA_anamiRTs_dr_fc[,c("NCD81", "NCD82", "NCD83")])/rowMeans(miRNA_anamiRTs_dr_fc[,c("ACD81", "ACD82", "ACD83")]))
  miRNA_bg_dr_fc$fc = log2(rowMeans(miRNA_bg_dr_fc[,c("NCD81", "NCD82", "NCD83")])/rowMeans(miRNA_bg_dr_fc[,c("ACD81", "ACD82", "ACD83")]))
  
  wilcox_test_inf <- wilcox.test(miRNA_inf_dr_fc$fc, miRNA_bg_dr_fc$fc)
  wilcox_test_ts <- wilcox.test(miRNA_ts_dr_fc$fc, miRNA_bg_dr_fc$fc)
  wilcox_test_tsnopct <- wilcox.test(miRNA_tsnopct_dr_fc$fc, miRNA_bg_dr_fc$fc)
  wilcox_test_both <- wilcox.test(miRNA_both_dr_fc$fc, miRNA_bg_dr_fc$fc)
  wilcox_test_bothnopct <- wilcox.test(miRNA_bothnopct_dr_fc$fc, miRNA_bg_dr_fc$fc)
  wilcox_test_miRDB <- wilcox.test(miRNA_miRDB_dr_fc$fc, miRNA_bg_dr_fc$fc)
  wilcox_test_miRDBTs <- wilcox.test(miRNA_miRDBTs_dr_fc$fc, miRNA_bg_dr_fc$fc)
  wilcox_test_miRWalk <- wilcox.test(miRNA_miRWalk_dr_fc$fc, miRNA_bg_dr_fc$fc)
  wilcox_test_miRWalkTs <- wilcox.test(miRNA_miRWalkTs_dr_fc$fc, miRNA_bg_dr_fc$fc)
  wilcox_test_anamiR <- wilcox.test(miRNA_anamiR_dr_fc$fc, miRNA_bg_dr_fc$fc)
  wilcox_test_anamiRTs <- wilcox.test(miRNA_anamiRTs_dr_fc$fc, miRNA_bg_dr_fc$fc)
  wilcox_test_bg <- wilcox.test(miRNA_bg_dr_fc$fc, miRNA_bg_dr_fc$fc)
  
  miRNA_inf_dr_fc$gene_name = rownames(miRNA_inf_dr_fc)
  miRNA_inf_dr_fc$type = 'Inferelator'
  miRNA_ts_dr_fc$gene_name = rownames(miRNA_ts_dr_fc)
  miRNA_ts_dr_fc$type = 'Targetscan conserved'
  miRNA_tsnopct_dr_fc$gene_name = rownames(miRNA_tsnopct_dr_fc)
  miRNA_tsnopct_dr_fc$type = 'Targetscan'
  miRNA_both_dr_fc$gene_name = rownames(miRNA_both_dr_fc)
  miRNA_both_dr_fc$type = 'Both conserved'
  miRNA_bothnopct_dr_fc$gene_name = rownames(miRNA_bothnopct_dr_fc)
  miRNA_bothnopct_dr_fc$type = 'Both'
  miRNA_miRDB_dr_fc$gene_name = rownames(miRNA_miRDB_dr_fc)
  miRNA_miRDB_dr_fc$type = 'miRDB'
  miRNA_miRDBTs_dr_fc$gene_name = rownames(miRNA_miRDBTs_dr_fc)
  miRNA_miRDBTs_dr_fc$type = 'miRDB & Targetscan'
  miRNA_miRWalk_dr_fc$gene_name = rownames(miRNA_miRWalk_dr_fc)
  miRNA_miRWalk_dr_fc$type = 'miRWalk'
  miRNA_miRWalkTs_dr_fc$gene_name = rownames(miRNA_miRWalkTs_dr_fc)
  miRNA_miRWalkTs_dr_fc$type = 'miRWalk & Targetscan'
  miRNA_anamiR_dr_fc$gene_name = rownames(miRNA_anamiR_dr_fc)
  miRNA_anamiR_dr_fc$type = 'anamiR'
  miRNA_anamiRTs_dr_fc$gene_name = rownames(miRNA_anamiRTs_dr_fc)
  miRNA_anamiRTs_dr_fc$type = 'anamiR & Targetscan'
  miRNA_bg_dr_fc$gene_name = rownames(miRNA_bg_dr_fc)
  miRNA_bg_dr_fc$type = 'Background'
  
  forplot = rbind(miRNA_ts_dr_fc, miRNA_tsnopct_dr_fc, miRNA_both_dr_fc, miRNA_bothnopct_dr_fc, #miRNA_inf_dr_fc, 
                  miRNA_miRDB_dr_fc, miRNA_miRDBTs_dr_fc, miRNA_miRWalk_dr_fc, miRNA_miRWalkTs_dr_fc, miRNA_anamiR_dr_fc, miRNA_anamiRTs_dr_fc,
                  miRNA_bg_dr_fc)
  forplot = forplot[,c('fc', 'type')]
  ggplot(forplot, aes(fc, colour = type)) +
    stat_ecdf(linewidth=1)+ 
    theme_classic() + #theme_classic
    coord_cartesian(xlim = c(-1, 1)) + # remove this line to get the nolimit plots
    #ylim(0, 1) + 
    xlab('Log Fold Change: NCD vs ACD') +
    ylab('Cumulative distribution') + 
    scale_color_manual(breaks = c('Background', 'Targetscan', 'Targetscan conserved', 'Both', 'Both conserved', 
                                  'miRDB', 'miRDB & Targetscan', 'miRWalk', 'miRWalk & Targetscan', 'anamiR', 'anamiR & Targetscan'),
                       labels = c(paste0("Background (n=", nrow(miRNA_bg_dr_fc), ")"),
                                  paste0("Targetscan (n=", nrow(miRNA_tsnopct_dr_fc), ", p=", formatC(wilcox_test_tsnopct$p.value, format = "e", digits = 2), ")"),
                                  paste0("Targetscan conserved (n=", nrow(miRNA_ts_dr_fc), ", p=", formatC(wilcox_test_ts$p.value, format = "e", digits = 2), ")"),
                                  #paste0("Inferelator (n=", nrow(miRNA_inf_dr_fc), ", p=", formatC(wilcox_test_inf$p.value, format = "e", digits = 2), ")"), 
                                  paste0("Inferelator (n=", nrow(miRNA_bothnopct_dr_fc), ", p=", formatC(wilcox_test_bothnopct$p.value, format = "e", digits = 2), ")"),
                                  paste0("Inferelator & TargetscanCon (n=", nrow(miRNA_both_dr_fc), ", p=", formatC(wilcox_test_both$p.value, format = "e", digits = 2), ")"),
                                  paste0("miRDB (n=", nrow(miRNA_miRDB_dr_fc), ", p=", formatC(wilcox_test_miRDB$p.value, format = "e", digits = 2), ")"),
                                  paste0("miRDB & TargetscanCon (n=", nrow(miRNA_miRDBTs_dr_fc), ", p=", formatC(wilcox_test_miRDBTs$p.value, format = "e", digits = 2), ")"),
                                  paste0("miRWalk (n=", nrow(miRNA_miRWalk_dr_fc), ", p=", formatC(wilcox_test_miRWalk$p.value, format = "e", digits = 2), ")"),
                                  paste0("miRWalk & TargetscanCon (n=", nrow(miRNA_miRWalkTs_dr_fc), ", p=", formatC(wilcox_test_miRWalkTs$p.value, format = "e", digits = 2), ")"), 
                                  paste0("anamiR (n=", nrow(miRNA_anamiR_dr_fc), ", p=", formatC(wilcox_test_anamiR$p.value, format = "e", digits = 2), ")"),
                                  paste0("anamiR & TargetscanCon (n=", nrow(miRNA_anamiRTs_dr_fc), ", p=", formatC(wilcox_test_anamiRTs$p.value, format = "e", digits = 2), ")")), 
                       values = c("black", "navyblue", "lightblue2", 'yellow4', 'yellow1', '#acc674', '#005c29', '#fa8072', '#800000', '#d7a1f9', '#7921b1'))+
    ggtitle(paste0(miRNA_name, " Targets\n", plot_title, " Fold Change"))#+
  #annotate("text", x=0.5, y=0.02, label=paste0('P-value = ', formatC(p_value_pos, format = "e", digits = 2)), size = 4)
  ggsave(
    paste0('/workdir/hz543/projects/Inferelator/mirna_gene_naiveOnly_noSP8_wt_prior/remake_figures/ct02_ms2_4ormore/benchmark/mir29/', 
                  "/cd8devproject_validatetargets_pct/", miRNA_name, "/", file_name, "_", file_path_sans_ext(basename(network_path)), "_lineplots_", miRNA_name, "_cd8devproject_ncdacd.pdf"),
    plot = last_plot(),
    device = "pdf",
    width = 7,
    height = 3.5,
    dpi = 300
  )
  
  
  png(file=paste0('make_figures/ct02_ms2_4ormore/benchmark/mir29/', 
                  "/cd8devproject_validatetargets_pct/", miRNA_name, "/", file_name, "_", file_path_sans_ext(basename(network_path)), "_lineplots_", miRNA_name, "_cd8devproject_ncdacd.png"), width=1000, height=1000, res=120)
  #par(mar=c(5.1, 4.1, 4.1, 18), xpd=TRUE)
  plot(as.list(environment(ecdf(miRNA_inf_dr_fc$fc))), pch = ".", col = "white", xlim=c(-1,1), 
       xlab = paste0("Log Fold Change: NCD vs ACD"), ylab = "Density", 
       cex.lab = 1.8, cex.axis = 1.8, cex.main = 2)
  #lines(as.list(environment(ecdf(miRNA_inf_dr_fc$fc))), col = '#cf1578', pch = ".", lwd = 3)
  lines(as.list(environment(ecdf(miRNA_ts_dr_fc$fc))), col = 'lightblue2', pch = ".", lwd = 3)
  lines(as.list(environment(ecdf(miRNA_tsnopct_dr_fc$fc))), col = 'navyblue', pch = ".", lwd = 3)
  lines(as.list(environment(ecdf(miRNA_miRDB_dr_fc$fc))), col = '#acc674', pch = ".", lwd = 3)
  lines(as.list(environment(ecdf(miRNA_miRDBTs_dr_fc$fc))), col = '#005c29', pch = ".", lwd = 3)
  lines(as.list(environment(ecdf(miRNA_miRWalk_dr_fc$fc))), col = '#fa8072', pch = ".", lwd = 3)
  lines(as.list(environment(ecdf(miRNA_miRWalkTs_dr_fc$fc))), col = '#800000', pch = ".", lwd = 3)
  lines(as.list(environment(ecdf(miRNA_anamiR_dr_fc$fc))), col = '#d7a1f9', pch = ".", lwd = 3)
  lines(as.list(environment(ecdf(miRNA_anamiRTs_dr_fc$fc))), col = '#7921b1', pch = ".", lwd = 3)
  lines(as.list(environment(ecdf(miRNA_both_dr_fc$fc))), col = 'yellow1', pch = ".", lwd = 3)
  lines(as.list(environment(ecdf(miRNA_bothnopct_dr_fc$fc))), col = 'yellow4', pch = ".", lwd = 3)
  lines(as.list(environment(ecdf(miRNA_bg_dr_fc$fc))), col = 'black', lwd = 3)
  abline(v=0, h=0.5, col = c("gray","gray"), lty = c(3,3))  
  legend(-1.05, 1.03, #inset=c(-0.2,0), 
         c(paste0("Background (n=", nrow(miRNA_bg_dr_fc), ")"),
           paste0("Targetscan (n=", nrow(miRNA_tsnopct_dr_fc), ", p=", formatC(wilcox_test_tsnopct$p.value, format = "e", digits = 2), ")"),
           paste0("Targetscan conserved (n=", nrow(miRNA_ts_dr_fc), ", p=", formatC(wilcox_test_ts$p.value, format = "e", digits = 2), ")"),
           #paste0("Inferelator (n=", nrow(miRNA_inf_dr_fc), ", p=", formatC(wilcox_test_inf$p.value, format = "e", digits = 2), ")"), 
           paste0("Both (n=", nrow(miRNA_bothnopct_dr_fc), ", p=", formatC(wilcox_test_bothnopct$p.value, format = "e", digits = 2), ")"),
           paste0("Both conserved (n=", nrow(miRNA_both_dr_fc), ", p=", formatC(wilcox_test_both$p.value, format = "e", digits = 2), ")"),
           paste0("miRDB (n=", nrow(miRNA_miRDB_dr_fc), ")"),
           paste0("miRDB&Targetscan (n=", nrow(miRNA_miRDBTs_dr_fc), ", p=", formatC(wilcox_test_tsnopct$p.value, format = "e", digits = 2), ")"),
           paste0("miRWalk (n=", nrow(miRNA_miRWalk_dr_fc), ", p=", formatC(wilcox_test_ts$p.value, format = "e", digits = 2), ")"),
           paste0("miRWalk&Targetscan (n=", nrow(miRNA_miRWalkTs_dr_fc), ", p=", formatC(wilcox_test_inf$p.value, format = "e", digits = 2), ")"), 
           paste0("anamiR (n=", nrow(miRNA_anamiR_dr_fc), ", p=", formatC(wilcox_test_bothnopct$p.value, format = "e", digits = 2), ")"),
           paste0("anamiR&Targetscan (n=", nrow(miRNA_anamiRTs_dr_fc), ", p=", formatC(wilcox_test_both$p.value, format = "e", digits = 2), ")")), 
         col = c("black", "navyblue", "lightblue2", 'yellow4', 'yellow1', '#acc674', '#005c29', '#fa8072', '#800000', '#d7a1f9', '#7921b1'),  # "#cf1578", 
         lty=c(1,1), 
         cex = 1.1, 
         lwd = 1.5,
         pt.cex = 10)
  title(main = paste0(miRNA_name, " Targets\n", plot_title, " Fold Change"), cex.main = 2)
  dev.off()
  
  
}




sample_info = read.table('inputs/mirna_exp_matrix/sample_info.txt', sep = '\t')

decay_rates_input = read.table('inputs/exonIntron_matrix/decay_rates.txt', header = T, row.names = 1)
cpm_wholegene_input = cd8dev_counts_wholegene_rel_cpm

#ms = 1
#shared_by_nets = 4

for (ms in c(2)){ #,2,3
  for (shared_by_nets in c(3)){ #,4
    
    bias = 5
    folderpath = paste0('analyze_network/ms', ms, '/bias', bias)
    
    for (seed_input in 'AGCACCA')
      for (tfaOpt_input in c("_TFmRNA")) # "", 
        for (iftsshared in c('')){ # , '_targetscansharedtargets'
          tryCatch({
            print(seed_input)
            
            plot_barplot_lineplots_cd8devproject(matrix = decay_rates_input, seed = seed_input, 
                                                 network_path = paste0(folderpath, '/sharedbyMorethan', shared_by_nets, 'netsFrom5seeds_pos/drcombat_ts_ct02_bias', bias, tfaOpt_input, '_cut01_sharedbyMorethan', shared_by_nets, 'nets_pos', iftsshared, '_sp.tsv'), 
                                                 plot_title = "Decay Rates", file_name = "dr")
            plot_barplot_lineplots_cd8devproject(matrix = cpm_wholegene_input, seed = seed_input, 
                                                 network_path = paste0(folderpath, '/sharedbyMorethan', shared_by_nets, 'netsFrom5seeds_pos/drcombat_ts_ct02_bias', bias, tfaOpt_input, '_cut01_sharedbyMorethan', shared_by_nets, 'nets_pos', iftsshared, '_sp.tsv'), 
                                                 plot_title = "Gene CPM", file_name = "gene")
          }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
        }
  }
}
