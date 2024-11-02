# get miRNA network
# consider only miRNA-gene regulations that appear four or five times in five runs; consider only positive regulations.

options(warn=1)

library(ggplot2)

for (ms in c(2)){
  
rm = ''
bias = 5
outputpath = paste0('analyze_network/ms', ms, '/bias', bias, '/prioritizeTargets')
dir.create(outputpath, showWarnings = FALSE, recursive = TRUE)

shared_by_nets = 0
num = shared_by_nets
#tfaOpt = c('', '_TFmRNA')
networkfolder = 'outputs'

tfaOpt = c('_TFmRNA')

sp_seed18 = read.table(paste0(networkfolder, '/mirnaSeed_tsonly_noSP8_wt_ct02_seed18', rm, '/bias', bias, tfaOpt, '_modelSize', ms, '/networks_targ0p05_SS50_bS5/Network0p05_', ms, 'tfsPerGene/tsmirnasexpressed1000_ct02_contextscore0to1_targGeneonly_bias', bias, tfaOpt, '_sp.tsv'), sep = '\t', header = T)
sp_seed26 = read.table(paste0(networkfolder, '/mirnaSeed_tsonly_noSP8_wt_ct02_seed26', rm, '/bias', bias, tfaOpt, '_modelSize', ms, '/networks_targ0p05_SS50_bS5/Network0p05_', ms, 'tfsPerGene/tsmirnasexpressed1000_ct02_contextscore0to1_targGeneonly_bias', bias, tfaOpt, '_sp.tsv'), sep = '\t', header = T)
sp_seed57 = read.table(paste0(networkfolder, '/mirnaSeed_tsonly_noSP8_wt_ct02_seed57', rm, '/bias', bias, tfaOpt, '_modelSize', ms, '/networks_targ0p05_SS50_bS5/Network0p05_', ms, 'tfsPerGene/tsmirnasexpressed1000_ct02_contextscore0to1_targGeneonly_bias', bias, tfaOpt, '_sp.tsv'), sep = '\t', header = T)
sp_seed7 = read.table(paste0(networkfolder, '/mirnaSeed_tsonly_noSP8_wt_ct02_seed7', rm, '/bias', bias, tfaOpt, '_modelSize', ms, '/networks_targ0p05_SS50_bS5/Network0p05_', ms, 'tfsPerGene/tsmirnasexpressed1000_ct02_contextscore0to1_targGeneonly_bias', bias, tfaOpt, '_sp.tsv'), sep = '\t', header = T)
sp_seed99 = read.table(paste0(networkfolder, '/mirnaSeed_tsonly_noSP8_wt_ct02_seed99', rm, '/bias', bias, tfaOpt, '_modelSize', ms, '/networks_targ0p05_SS50_bS5/Network0p05_', ms, 'tfsPerGene/tsmirnasexpressed1000_ct02_contextscore0to1_targGeneonly_bias', bias, tfaOpt, '_sp.tsv'), sep = '\t', header = T)

sp_seed18$direction = NA
for (i in 1:nrow(sp_seed18)){
  if (sp_seed18[i,'SignedQuantile'] > 0){
    sp_seed18[i,'direction'] = 'pos'
  } else if (sp_seed18[i,'SignedQuantile'] < 0){
    sp_seed18[i,'direction'] = 'neg'
  }
}
sp_seed18$'TF->Target' = paste(sp_seed18$TF, sp_seed18$Target, sep = '->')
sp_seed18$'TF->Target,direction' = paste(sp_seed18$'TF->Target', sp_seed18$direction, sep = ',')
sp_seed18$TF = NULL
sp_seed18$Target = NULL
sp_seed18$stroke.width = NULL
sp_seed18$stroke.dasharray = NULL
sp_seed18$stroke = NULL
sp_seed18$`TF->Target` = NULL
sp_seed18_pos = sp_seed18[sp_seed18$direction == 'pos', ]
sp_seed18_neg = sp_seed18[sp_seed18$direction == 'neg', ]


sp_seed7$direction = NA
for (i in 1:nrow(sp_seed7)){
  if (sp_seed7[i,'SignedQuantile'] > 0){
    sp_seed7[i,'direction'] = 'pos'
  } else if (sp_seed7[i,'SignedQuantile'] < 0){
    sp_seed7[i,'direction'] = 'neg'
  }
}
sp_seed7$'TF->Target' = paste(sp_seed7$TF, sp_seed7$Target, sep = '->')
sp_seed7$'TF->Target,direction' = paste(sp_seed7$'TF->Target', sp_seed7$direction, sep = ',')
sp_seed7$TF = NULL
sp_seed7$Target = NULL
sp_seed7$stroke.width = NULL
sp_seed7$stroke.dasharray = NULL
sp_seed7$stroke = NULL
sp_seed7$`TF->Target` = NULL
sp_seed7_pos = sp_seed7[sp_seed7$direction == 'pos', ]
sp_seed7_neg = sp_seed7[sp_seed7$direction == 'neg', ]

sp_seed99$direction = NA
for (i in 1:nrow(sp_seed99)){
  if (sp_seed99[i,'SignedQuantile'] > 0){
    sp_seed99[i,'direction'] = 'pos'
  } else if (sp_seed99[i,'SignedQuantile'] < 0){
    sp_seed99[i,'direction'] = 'neg'
  }
}
sp_seed99$'TF->Target' = paste(sp_seed99$TF, sp_seed99$Target, sep = '->')
sp_seed99$'TF->Target,direction' = paste(sp_seed99$'TF->Target', sp_seed99$direction, sep = ',')
sp_seed99$TF = NULL
sp_seed99$Target = NULL
sp_seed99$stroke.width = NULL
sp_seed99$stroke.dasharray = NULL
sp_seed99$stroke = NULL
sp_seed99$`TF->Target` = NULL
sp_seed99_pos = sp_seed99[sp_seed99$direction == 'pos', ]
sp_seed99_neg = sp_seed99[sp_seed99$direction == 'neg', ]

sp_seed26$direction = NA
for (i in 1:nrow(sp_seed26)){
  if (sp_seed26[i,'SignedQuantile'] > 0){
    sp_seed26[i,'direction'] = 'pos'
  } else if (sp_seed26[i,'SignedQuantile'] < 0){
    sp_seed26[i,'direction'] = 'neg'
  }
}
sp_seed26$'TF->Target' = paste(sp_seed26$TF, sp_seed26$Target, sep = '->')
sp_seed26$'TF->Target,direction' = paste(sp_seed26$'TF->Target', sp_seed26$direction, sep = ',')
sp_seed26$TF = NULL
sp_seed26$Target = NULL
sp_seed26$stroke.width = NULL
sp_seed26$stroke.dasharray = NULL
sp_seed26$stroke = NULL
sp_seed26$`TF->Target` = NULL
sp_seed26_pos = sp_seed26[sp_seed26$direction == 'pos', ]
sp_seed26_neg = sp_seed26[sp_seed26$direction == 'neg', ]

sp_seed57$direction = NA
for (i in 1:nrow(sp_seed57)){
  if (sp_seed57[i,'SignedQuantile'] > 0){
    sp_seed57[i,'direction'] = 'pos'
  } else if (sp_seed57[i,'SignedQuantile'] < 0){
    sp_seed57[i,'direction'] = 'neg'
  }
}
sp_seed57$'TF->Target' = paste(sp_seed57$TF, sp_seed57$Target, sep = '->')
sp_seed57$'TF->Target,direction' = paste(sp_seed57$'TF->Target', sp_seed57$direction, sep = ',')
sp_seed57$TF = NULL
sp_seed57$Target = NULL
sp_seed57$stroke.width = NULL
sp_seed57$stroke.dasharray = NULL
sp_seed57$stroke = NULL
sp_seed57$`TF->Target` = NULL
sp_seed57_pos = sp_seed57[sp_seed57$direction == 'pos', ]
sp_seed57_neg = sp_seed57[sp_seed57$direction == 'neg', ]


res_table = table(c(sp_seed18$`TF->Target,direction`, sp_seed7$`TF->Target,direction`, sp_seed99$`TF->Target,direction`, sp_seed26$`TF->Target,direction`, sp_seed57$`TF->Target,direction`))
res_table_pos = table(c(sp_seed18_pos$`TF->Target,direction`, sp_seed7_pos$`TF->Target,direction`, sp_seed99_pos$`TF->Target,direction`, sp_seed26_pos$`TF->Target,direction`, sp_seed57_pos$`TF->Target,direction`))
res_table_neg = table(c(sp_seed18_neg$`TF->Target,direction`, sp_seed7_neg$`TF->Target,direction`, sp_seed99_neg$`TF->Target,direction`, sp_seed26_neg$`TF->Target,direction`, sp_seed57_neg$`TF->Target,direction`))
#sp_all = merge(sp_seed18, sp_seed7, sp_seed99, sp_seed26, sp_seed57, by = 'TF->Target,direction', all = T)

MyMerge <- function(x, y){
  df <- merge(x, y, by= "TF->Target,direction", all = TRUE)
  return(df)
}
sp_all <- Reduce(MyMerge, list(sp_seed18, sp_seed7, sp_seed99, sp_seed26, sp_seed57))
sp_all_pos <- Reduce(MyMerge, list(sp_seed18_pos, sp_seed7_pos, sp_seed99_pos, sp_seed26_pos, sp_seed57_pos))
sp_all_neg <- Reduce(MyMerge, list(sp_seed18_neg, sp_seed7_neg, sp_seed99_neg, sp_seed26_neg, sp_seed57_neg))

colnames(sp_all_pos) = c("TF->Target,direction", 
                         "netTFmRNA1_SignedQuantile", "netTFmRNA1_NonzeroSubsamples.Stability", "netTFmRNA1_pCorr", "netTFmRNA1_direction", 
                         "netTFmRNA2_SignedQuantile", "netTFmRNA2_NonzeroSubsamples.Stability", "netTFmRNA2_pCorr", "netTFmRNA2_direction", 
                         "netTFmRNA3_SignedQuantile", "netTFmRNA3_NonzeroSubsamples.Stability", "netTFmRNA3_pCorr", "netTFmRNA3_direction", 
                         "netTFmRNA4_SignedQuantile", "netTFmRNA4_NonzeroSubsamples.Stability", "netTFmRNA4_pCorr", "netTFmRNA4_direction", 
                         "netTFmRNA5_SignedQuantile", "netTFmRNA5_NonzeroSubsamples.Stability", "netTFmRNA5_pCorr", "netTFmRNA5_direction")
sp_all_pos$netTFmRNA1_direction = NULL
sp_all_pos$netTFmRNA2_direction = NULL
sp_all_pos$netTFmRNA3_direction = NULL
sp_all_pos$netTFmRNA4_direction = NULL
sp_all_pos$netTFmRNA5_direction = NULL

sp_all_pos$freqTFmRNA = NA
for (i in 1:length(res_table_pos)){
  sp_all_pos[sp_all_pos$`TF->Target,direction` == names(res_table_pos[i]), 'freqTFmRNA'] = res_table_pos[i]
}

pos_regs_TFmRNA = sp_all_pos


# from Targetscan, find miRNAs' seed sequence (used to find targets) 
mir_family_info <- read.table("targetscan/mouse/miR_Family_Info.txt", head = T, fill = T, sep = "\t") 
hsa_mir_family_info <- mir_family_info[substr(mir_family_info[,'MiRBase.ID'] ,1 ,3) == "mmu", ] 
hsa_con_mir_family_info <- hsa_mir_family_info[hsa_mir_family_info[,'Family.Conservation.'] >= 2, ] 
seed_con <- hsa_con_mir_family_info[,c('Seed.m8', 'MiRBase.ID')] 
mirna_seeds_con <- aggregate(seed_con[,2], by = list(seed = seed_con$Seed), FUN = paste) 
colnames(mirna_seeds_con) <- c("seed", "miRNA") 
mirna_seeds_con$miRNA_name = c('miR-106ab,17,20ab,93-5p,6383','miR-141,200a-3p','miR-132,212-3p','miR-451a','miR-490-3p','miR-191-5p','miR-124-3p.1','miR-18ab-5p','miR-291a,294,295,302abd-3p','miR-200bc,429-3p','miR-216a-5p','miR-216b-5p','miR-365-3p','miR-101a-3p.1','miR-144-3p','miR-181abcd-5p','miR-140-3p.2,miR-497b','miR-100,99ab-5p','miR-10ab-5p','miR-217-5p','miR-193ab-3p','miR-383-5p.2','miR-29abc-3p','miR-15ab,16,1907,195a,322,497a-5p,miR-195b,6342,6353,6419','miR-129-1,2-3p','miR-22-3p','miR-21a,590-5p,miR-21c','miR-196ab-5p','miR-130ab,301ab-3p,miR-130c,6341,6389,721','miR-302c-3p','miR-140-5p,miR-876-3p','miR-142a-5p','miR-425-5p,miR-489-3p','miR-183-5p','miR-135ab-5p','miR-455-5p,miR-5129-3p','miR-25,363,367,92ab-3p,miR-32-5p','miR-128-3p,miR-6539','miR-802-5p','miR-199ab-3p','miR-455-3p.1','miR-148ab,152-3p','miR-140-3p.1','miR-338-3p','miR-199ab-5p','miR-125ab,351-5p,miR-6367,6394','miR-205-5p','miR-212-5p','miR-551b-3p','miR-126a-3p.1','miR-187-3p','miR-139-5p','miR-150-5p,miR-5127','miR-9-5p','miR-203-3p.1','miR-146ab-5p','miR-143-3p','let-7abcdefgik,98-5p,miR-1961','miR-190ab-5p','miR-383-5p.1','miR-219a-5p','miR-103,107-3p','miR-214-5p','miR-221,222-3p,miR-1928','miR-138-5p','miR-7ab-5p','miR-1a,206-3p,miR-1957b,6349,6382','miR-184-3p','miR-122-5p','miR-31-5p','miR-34abc,449ac-5p,-miR-449b','miR-24-3p,miR-5124b,6361,6369,6410,6413','miR-193a-5p','miR-30abcde,384-5p','miR-194-5p','miR-126a-3p.2','miR-142a-3p.1','miR-223-3p','miR-19ab-3p','miR-208ab-3p','miR-499-5p','miR-124,5624-3p,miR-6540-5p','miR-155-5p','miR-101ab-3p','miR-142a-3p.2','miR-137-3p','miR-26ab-5p','miR-27ab-3p','miR-23ab-3p','miR-145a-5p,miR-145b','miR-204,211-5p,miR-7670-3p','miR-202-5p','miR-203-3p.2','miR-192,215-5p','miR-455-3p.2,miR-682','miR-153-3p','miR-33-5p','miR-183-5p.2','miR-133a-3p.1','miR-147-3p','miR-210-3p','miR-218-5p,miR-7002-3p','miR-182-5p','miR-96-5p','miR-133ab-3p,miR-133c','miR-375-3p','miR-129-5p')

pos_regs_TFmRNA$miRNA_seed = gsub('->.*', '', pos_regs_TFmRNA$`TF->Target,direction`)
pos_regs_TFmRNA$target = gsub('.*->', '', pos_regs_TFmRNA$`TF->Target,direction`)
pos_regs_TFmRNA$target = gsub(',.*', '', pos_regs_TFmRNA$target)
pos_regs_TFmRNA$miRNA_name = NA
for(i in 1:nrow(pos_regs_TFmRNA)){
  pos_regs_TFmRNA[i, 'miRNA_name'] = mirna_seeds_con[mirna_seeds_con$seed == pos_regs_TFmRNA[i, 'miRNA_seed'],'miRNA_name']
}
pos_regs_TFmRNA = pos_regs_TFmRNA[, c("TF->Target,direction", "miRNA_seed", "target", "miRNA_name", "freqTFmRNA",
                                "netTFmRNA1_SignedQuantile", "netTFmRNA1_NonzeroSubsamples.Stability", "netTFmRNA1_pCorr", 
                                "netTFmRNA2_SignedQuantile", "netTFmRNA2_NonzeroSubsamples.Stability", "netTFmRNA2_pCorr", 
                                "netTFmRNA3_SignedQuantile", "netTFmRNA3_NonzeroSubsamples.Stability", "netTFmRNA3_pCorr", 
                                "netTFmRNA4_SignedQuantile", "netTFmRNA4_NonzeroSubsamples.Stability", "netTFmRNA4_pCorr", 
                                "netTFmRNA5_SignedQuantile", "netTFmRNA5_NonzeroSubsamples.Stability", "netTFmRNA5_pCorr")]




# TFmRNA option more than 4

pos_regs_TFmRNA_freqTFmRNAmorethan4 = pos_regs_TFmRNA[pos_regs_TFmRNA$freqTFmRNA >= 4, ]
pos_regs_TFmRNA_freqTFmRNAmorethan4 = pos_regs_TFmRNA_freqTFmRNAmorethan4[!is.na(pos_regs_TFmRNA_freqTFmRNAmorethan4$miRNA_seed),]

ts_nopct = read.table('inputs/prior/ts_ct02_contextscore.tsv')
ts = read.table('inputs/prior/ts_ct02_pct75_b.tsv')

pos_regs_TFmRNA_freqTFmRNAmorethan4$ts_nopct = NA
for(i in 1:nrow(pos_regs_TFmRNA_freqTFmRNAmorethan4)){
  if (!pos_regs_TFmRNA_freqTFmRNAmorethan4[i,'target'] %in% rownames(ts_nopct)){
    pos_regs_TFmRNA_freqTFmRNAmorethan4[i, 'ts_nopct'] = FALSE
  } else if (ts_nopct[pos_regs_TFmRNA_freqTFmRNAmorethan4[i,'target'], pos_regs_TFmRNA_freqTFmRNAmorethan4[i,'miRNA_seed']] != 0){
    pos_regs_TFmRNA_freqTFmRNAmorethan4[i, 'ts_nopct'] = TRUE
  } else {
    pos_regs_TFmRNA_freqTFmRNAmorethan4[i, 'ts_nopct'] = FALSE
  }
}

pos_regs_TFmRNA_freqTFmRNAmorethan4$ts = NA
for(i in 1:nrow(pos_regs_TFmRNA_freqTFmRNAmorethan4)){
  if (pos_regs_TFmRNA_freqTFmRNAmorethan4[i,'target'] %in% rownames(ts)){
    if (ts[pos_regs_TFmRNA_freqTFmRNAmorethan4[i,'target'], pos_regs_TFmRNA_freqTFmRNAmorethan4[i,'miRNA_seed']] != 0){
      pos_regs_TFmRNA_freqTFmRNAmorethan4[i, 'ts'] = TRUE
    } else {
      pos_regs_TFmRNA_freqTFmRNAmorethan4[i, 'ts'] = FALSE
    }
  }
}
pos_regs_TFmRNA_freqTFmRNAmorethan4 = pos_regs_TFmRNA_freqTFmRNAmorethan4[, c("TF->Target,direction", "miRNA_seed", "target", "miRNA_name", "freqTFmRNA", "ts_nopct", "ts",
                                                "netTFmRNA1_SignedQuantile", "netTFmRNA1_NonzeroSubsamples.Stability", "netTFmRNA1_pCorr", 
                                                "netTFmRNA2_SignedQuantile", "netTFmRNA2_NonzeroSubsamples.Stability", "netTFmRNA2_pCorr", 
                                                "netTFmRNA3_SignedQuantile", "netTFmRNA3_NonzeroSubsamples.Stability", "netTFmRNA3_pCorr", 
                                                "netTFmRNA4_SignedQuantile", "netTFmRNA4_NonzeroSubsamples.Stability", "netTFmRNA4_pCorr", 
                                                "netTFmRNA5_SignedQuantile", "netTFmRNA5_NonzeroSubsamples.Stability", "netTFmRNA5_pCorr")]

write.table(pos_regs_TFmRNA_freqTFmRNAmorethan4, paste0(outputpath, '/pos_regs_freqTFmRNAmorethan4.tsv'),
            sep = '\t', quote = F, row.names = F)  


}


