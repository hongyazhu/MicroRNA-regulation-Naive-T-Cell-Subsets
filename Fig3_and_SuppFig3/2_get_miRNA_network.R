# get network shared by multiple runs - try sharing by more than 1, 2, 3, 4 runs

outputfolder = 'analyze_network'
networkfolder = 'outputs'

getsharednets <- function(ms, bias, tfaOpt, rm){
  
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
  sp_seed57_pos = sp_seed57[sp_seed57$direction == 'pos', ]
  sp_seed57_neg = sp_seed57[sp_seed57$direction == 'neg', ]
  
  #table(sort(table(c(sp_seed18$`TF->Target,direction`, sp_seed7$`TF->Target,direction`, sp_seed99$`TF->Target,direction`, sp_seed26$`TF->Target,direction`, sp_seed57$`TF->Target,direction`)), decreasing = T) == 5)
  print('Pos + neg: shared by more than 4 nets:')
  table(sort(table(c(sp_seed18$`TF->Target,direction`, sp_seed7$`TF->Target,direction`, sp_seed99$`TF->Target,direction`, sp_seed26$`TF->Target,direction`, sp_seed57$`TF->Target,direction`)), decreasing = T) > 4)
  print('Pos + neg: shared by more than 3 nets:')
  table(sort(table(c(sp_seed18$`TF->Target,direction`, sp_seed7$`TF->Target,direction`, sp_seed99$`TF->Target,direction`, sp_seed26$`TF->Target,direction`, sp_seed57$`TF->Target,direction`)), decreasing = T) > 3)
  print('Pos + neg: shared by more than 2 nets:')
  table(sort(table(c(sp_seed18$`TF->Target,direction`, sp_seed7$`TF->Target,direction`, sp_seed99$`TF->Target,direction`, sp_seed26$`TF->Target,direction`, sp_seed57$`TF->Target,direction`)), decreasing = T) > 2)
  print('Pos + neg: shared by more than 1 nets:')
  table(sort(table(c(sp_seed18$`TF->Target,direction`, sp_seed7$`TF->Target,direction`, sp_seed99$`TF->Target,direction`, sp_seed26$`TF->Target,direction`, sp_seed57$`TF->Target,direction`)), decreasing = T) > 1)
  print('Pos + neg: shared by more than 0 nets:')
  table(sort(table(c(sp_seed18$`TF->Target,direction`, sp_seed7$`TF->Target,direction`, sp_seed99$`TF->Target,direction`, sp_seed26$`TF->Target,direction`, sp_seed57$`TF->Target,direction`)), decreasing = T) > 0)
  
  
  print('Pos: shared by more than 4 nets:')
  table(sort(table(c(sp_seed18_pos$`TF->Target,direction`, sp_seed7_pos$`TF->Target,direction`, sp_seed99_pos$`TF->Target,direction`, sp_seed26_pos$`TF->Target,direction`, sp_seed57_pos$`TF->Target,direction`)), decreasing = T) > 4)
  print('Pos: shared by more than 3 nets:')
  table(sort(table(c(sp_seed18_pos$`TF->Target,direction`, sp_seed7_pos$`TF->Target,direction`, sp_seed99_pos$`TF->Target,direction`, sp_seed26_pos$`TF->Target,direction`, sp_seed57_pos$`TF->Target,direction`)), decreasing = T) > 3)
  print('Pos: shared by more than 2 nets:')
  table(sort(table(c(sp_seed18_pos$`TF->Target,direction`, sp_seed7_pos$`TF->Target,direction`, sp_seed99_pos$`TF->Target,direction`, sp_seed26_pos$`TF->Target,direction`, sp_seed57_pos$`TF->Target,direction`)), decreasing = T) > 2)
  print('Pos: shared by more than 1 nets:')
  table(sort(table(c(sp_seed18_pos$`TF->Target,direction`, sp_seed7_pos$`TF->Target,direction`, sp_seed99_pos$`TF->Target,direction`, sp_seed26_pos$`TF->Target,direction`, sp_seed57_pos$`TF->Target,direction`)), decreasing = T) > 1)
  print('Pos: shared by more than 0 nets:')
  table(sort(table(c(sp_seed18_pos$`TF->Target,direction`, sp_seed7_pos$`TF->Target,direction`, sp_seed99_pos$`TF->Target,direction`, sp_seed26_pos$`TF->Target,direction`, sp_seed57_pos$`TF->Target,direction`)), decreasing = T) > 0)
  
  
  print('Neg: shared by more than 4 nets:')
  table(sort(table(c(sp_seed18_neg$`TF->Target,direction`, sp_seed7_neg$`TF->Target,direction`, sp_seed99_neg$`TF->Target,direction`, sp_seed26_neg$`TF->Target,direction`, sp_seed57_neg$`TF->Target,direction`)), decreasing = T) > 4)
  print('Neg: shared by more than 3 nets:')
  table(sort(table(c(sp_seed18_neg$`TF->Target,direction`, sp_seed7_neg$`TF->Target,direction`, sp_seed99_neg$`TF->Target,direction`, sp_seed26_neg$`TF->Target,direction`, sp_seed57_neg$`TF->Target,direction`)), decreasing = T) > 3)
  print('Neg: shared by more than 2 nets:')
  table(sort(table(c(sp_seed18_neg$`TF->Target,direction`, sp_seed7_neg$`TF->Target,direction`, sp_seed99_neg$`TF->Target,direction`, sp_seed26_neg$`TF->Target,direction`, sp_seed57_neg$`TF->Target,direction`)), decreasing = T) > 2)
  print('Neg: shared by more than 1 nets:')
  table(sort(table(c(sp_seed18_neg$`TF->Target,direction`, sp_seed7_neg$`TF->Target,direction`, sp_seed99_neg$`TF->Target,direction`, sp_seed26_neg$`TF->Target,direction`, sp_seed57_neg$`TF->Target,direction`)), decreasing = T) > 1)
  print('Neg: shared by more than 0 nets:')
  table(sort(table(c(sp_seed18_neg$`TF->Target,direction`, sp_seed7_neg$`TF->Target,direction`, sp_seed99_neg$`TF->Target,direction`, sp_seed26_neg$`TF->Target,direction`, sp_seed57_neg$`TF->Target,direction`)), decreasing = T) > 0)
  
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
  
  #num = 4
  for (num in 1:4){
    dir.create(paste0(outputfolder, '/ms', ms, '/bias', bias, '/sharedbyMorethan', num, 'netsFrom5seeds/'), recursive = TRUE, showWarnings = FALSE)
    dir.create(paste0(outputfolder, '/ms', ms, '/bias', bias, '/sharedbyMorethan', num, 'netsFrom5seeds_pos/'), recursive = TRUE, showWarnings = FALSE)
    dir.create(paste0(outputfolder, '/ms', ms, '/bias', bias, '/sharedbyMorethan', num, 'netsFrom5seeds_neg/'), recursive = TRUE, showWarnings = FALSE)
  }
  
  for (num in 1:4){
    
    sp_interest_names = names(res_table[res_table > num])
    
    sp_interest = sp_all[sp_all$`TF->Target,direction` %in% sp_interest_names, ]
    
    sp_interest_less = sp_interest[, c('TF->Target,direction', 'TF.x')]
    sp_interest_less$TF = gsub('->.*', '', sp_interest_less$`TF->Target,direction`)
    sp_interest_less$Target = gsub(',.*', '', gsub('.*->', '', sp_interest_less$`TF->Target,direction`))
    sp_interest_less$direction = gsub('.*,', '', sp_interest_less$`TF->Target,direction`)
    sp_interest_less$SignedQuantile = NA
    for (i in 1:nrow(sp_interest_less)){
      if (sp_interest_less[i,'direction'] == 'pos'){
        sp_interest_less[i,'SignedQuantile'] = 1
      } else if (sp_interest_less[i,'direction'] == 'neg'){
        sp_interest_less[i,'SignedQuantile'] = -1
      }
    }
    sp_interest_format = sp_interest_less[, c('TF', 'Target', 'SignedQuantile')]
    print(num)
    print(nrow(sp_interest_format))
    write.table(sp_interest_format, paste0(outputfolder, '/ms', ms, '/bias', bias, '/sharedbyMorethan', num, 'netsFrom5seeds/drcombat_ts_ct02_bias', bias, tfaOpt, rm, '_cut01_sharedbyMorethan', num, 'nets_sp.tsv'), sep = '\t', quote = F, row.names = F)  
    
    
    notsp = data.frame(matrix(0, ncol = length(unique(sp_interest_format$TF)), nrow = length(unique(sp_interest_format$Target))))
    colnames(notsp) = unique(sp_interest_format$TF)
    rownames(notsp) = unique(sp_interest_format$Target)
    
    for (i in 1:nrow(sp_interest_format))
      if (sp_interest_format[i, 'SignedQuantile'] > 0){
        notsp[sp_interest_format[i,'Target'], sp_interest_format[i,'TF']] = 1
      } else if (sp_interest_format[i, 'SignedQuantile'] < 0){
        notsp[sp_interest_format[i,'Target'], sp_interest_format[i,'TF']] = -1
      }
    
    notsp_new <- notsp[ order(row.names(notsp)), ]
    notsp_new <- notsp[ , order(colnames(notsp))]
    
    
    write.table(notsp_new, paste0(outputfolder, '/ms', ms, '/bias', bias, '/sharedbyMorethan', num, 'netsFrom5seeds/drcombat_ts_ct02_bias', bias, tfaOpt, rm, '_cut01_sharedbyMorethan', num, 'nets.tsv'), sep = '\t', quote = F, col.names=NA)
  }
  
  
  
  for (num in 1:4){
    
    sp_interest_names = names(res_table_pos[res_table_pos > num])
    
    sp_interest = sp_all_pos[sp_all_pos$`TF->Target,direction` %in% sp_interest_names, ]
    
    sp_interest_less = sp_interest[, c('TF->Target,direction', 'TF.x')]
    sp_interest_less$TF = gsub('->.*', '', sp_interest_less$`TF->Target,direction`)
    sp_interest_less$Target = gsub(',.*', '', gsub('.*->', '', sp_interest_less$`TF->Target,direction`))
    sp_interest_less$direction = gsub('.*,', '', sp_interest_less$`TF->Target,direction`)
    sp_interest_less$SignedQuantile = NA
    for (i in 1:nrow(sp_interest_less)){
      if (sp_interest_less[i,'direction'] == 'pos'){
        sp_interest_less[i,'SignedQuantile'] = 1
      } else if (sp_interest_less[i,'direction'] == 'neg'){
        sp_interest_less[i,'SignedQuantile'] = -1
      }
    }
    sp_interest_format = sp_interest_less[, c('TF', 'Target', 'SignedQuantile')]
    print(num)
    print(nrow(sp_interest_format))
    write.table(sp_interest_format, paste0(outputfolder, '/ms', ms, '/bias', bias, '/sharedbyMorethan', num, 'netsFrom5seeds_pos/drcombat_ts_ct02_bias', bias, tfaOpt, rm, '_cut01_sharedbyMorethan', num, 'nets_pos_sp.tsv'), sep = '\t', quote = F, row.names = F)  
    
    
    notsp = data.frame(matrix(0, ncol = length(unique(sp_interest_format$TF)), nrow = length(unique(sp_interest_format$Target))))
    colnames(notsp) = unique(sp_interest_format$TF)
    rownames(notsp) = unique(sp_interest_format$Target)
    
    for (i in 1:nrow(sp_interest_format))
      if (sp_interest_format[i, 'SignedQuantile'] > 0){
        notsp[sp_interest_format[i,'Target'], sp_interest_format[i,'TF']] = 1
      } else if (sp_interest_format[i, 'SignedQuantile'] < 0){
        notsp[sp_interest_format[i,'Target'], sp_interest_format[i,'TF']] = -1
      }
    
    notsp_new <- notsp[ order(row.names(notsp)), ]
    notsp_new <- notsp[ , order(colnames(notsp))]
    
    
    write.table(notsp_new, paste0(outputfolder, '/ms', ms, '/bias', bias, '/sharedbyMorethan', num, 'netsFrom5seeds_pos/drcombat_ts_ct02_bias', bias, tfaOpt, rm, '_cut01_sharedbyMorethan', num, 'nets_pos.tsv'), sep = '\t', quote = F, col.names=NA)
  }
  
  
  for (num in 1:4){
    
    sp_interest_names = names(res_table_neg[res_table_neg > num])
    
    if (length(sp_interest_names) != 0){
      sp_interest = sp_all_neg[sp_all_neg$`TF->Target,direction` %in% sp_interest_names, ]
      
      sp_interest_less = sp_interest[, c('TF->Target,direction', 'TF.x')]
      sp_interest_less$TF = gsub('->.*', '', sp_interest_less$`TF->Target,direction`)
      sp_interest_less$Target = gsub(',.*', '', gsub('.*->', '', sp_interest_less$`TF->Target,direction`))
      sp_interest_less$direction = gsub('.*,', '', sp_interest_less$`TF->Target,direction`)
      sp_interest_less$SignedQuantile = NA
      for (i in 1:nrow(sp_interest_less)){
        if (sp_interest_less[i,'direction'] == 'pos'){
          sp_interest_less[i,'SignedQuantile'] = 1
        } else if (sp_interest_less[i,'direction'] == 'neg'){
          sp_interest_less[i,'SignedQuantile'] = -1
        }
      }
      sp_interest_format = sp_interest_less[, c('TF', 'Target', 'SignedQuantile')]
      print(num)
      print(nrow(sp_interest_format))
      write.table(sp_interest_format, paste0(outputfolder, '/ms', ms, '/bias', bias, '/sharedbyMorethan', num, 'netsFrom5seeds_neg/drcombat_ts_ct02_bias', bias, tfaOpt, rm, '_cut01_sharedbyMorethan', num, 'nets_neg_sp.tsv'), sep = '\t', quote = F, row.names = F)  
      
      
      notsp = data.frame(matrix(0, ncol = length(unique(sp_interest_format$TF)), nrow = length(unique(sp_interest_format$Target))))
      colnames(notsp) = unique(sp_interest_format$TF)
      rownames(notsp) = unique(sp_interest_format$Target)
      
      for (i in 1:nrow(sp_interest_format))
        if (sp_interest_format[i, 'SignedQuantile'] > 0){
          notsp[sp_interest_format[i,'Target'], sp_interest_format[i,'TF']] = 1
        } else if (sp_interest_format[i, 'SignedQuantile'] < 0){
          notsp[sp_interest_format[i,'Target'], sp_interest_format[i,'TF']] = -1
        }
      
      notsp_new <- notsp[ order(row.names(notsp)), ]
      notsp_new <- notsp[ , order(colnames(notsp))]
      
      
      write.table(notsp_new, paste0(outputfolder, '/ms', ms, '/bias', bias, '/sharedbyMorethan', num, 'netsFrom5seeds_neg/drcombat_ts_ct02_bias', bias, tfaOpt, rm, '_cut01_sharedbyMorethan', num, 'nets_neg.tsv'), sep = '\t', quote = F, col.names=NA)
    } else {
      print(paste0('sharedby', num, ': no neg regs'))
    }
  }
  
  
  ### shared with targetscan (for shared nets pos only)
  
  ts_ct02_pct75_b = read.table('/workdir/hz543/projects/Inferelator/mirna_gene_naiveOnly_noSP8_wt_prior/inputs/prior/ts_ct02_pct75_b.tsv')
  #num = 2
  
  for (num in 1:4){
    sp_interest_names = names(res_table_pos[res_table_pos > num])
    
    sp_interest = sp_all_pos[sp_all_pos$`TF->Target,direction` %in% sp_interest_names, ]
    
    sp_interest_less = sp_interest[, c('TF->Target,direction', 'TF.x')]
    sp_interest_less$TF = gsub('->.*', '', sp_interest_less$`TF->Target,direction`)
    sp_interest_less$Target = gsub(',.*', '', gsub('.*->', '', sp_interest_less$`TF->Target,direction`))
    sp_interest_less$direction = gsub('.*,', '', sp_interest_less$`TF->Target,direction`)
    sp_interest_less$SignedQuantile = NA
    for (i in 1:nrow(sp_interest_less)){
      if (sp_interest_less[i,'direction'] == 'pos'){
        sp_interest_less[i,'SignedQuantile'] = 1
      } else if (sp_interest_less[i,'direction'] == 'neg'){
        sp_interest_less[i,'SignedQuantile'] = -1
      }
    }
    sp_interest_format = sp_interest_less[, c('TF', 'Target', 'SignedQuantile')]
    print(num)
    print(nrow(sp_interest_format))
    
    print(table(unique(sp_interest_format$TF) %in% colnames(ts_ct02_pct75_b))) # all true - all TFs(miRNAs) are listed in targetscan
    print(table(unique(sp_interest_format$Target) %in% rownames(ts_ct02_pct75_b))) # not all true - not all targets are listed in targetscan
    
    sp_interest_format_targetscan = sp_interest_format[sp_interest_format$Target %in% rownames(ts_ct02_pct75_b),]
    print(table(unique(sp_interest_format_targetscan$Target) %in% rownames(ts_ct02_pct75_b))) # all true 
    
    shared_targets = c()
    for (i in 1:nrow(sp_interest_format_targetscan)){
      if (ts_ct02_pct75_b[sp_interest_format_targetscan[i,'Target'], sp_interest_format_targetscan[i,'TF']] != 0){
        shared_targets = c(shared_targets, i)
      }
    }
    sp_interest_format_targetscansharedtargets = sp_interest_format_targetscan[shared_targets, ]
    write.table(sp_interest_format_targetscansharedtargets, paste0(outputfolder, '/ms', ms, '/bias', bias, '/sharedbyMorethan', num, 'netsFrom5seeds_pos/drcombat_ts_ct02_bias', bias, tfaOpt, rm, '_cut01_sharedbyMorethan', num, 'nets_pos_targetscansharedtargets_sp.tsv'), sep = '\t', quote = F, row.names = F)  
    
    
    notsp = data.frame(matrix(0, ncol = length(unique(sp_interest_format_targetscansharedtargets$TF)), nrow = length(unique(sp_interest_format_targetscansharedtargets$Target))))
    colnames(notsp) = unique(sp_interest_format_targetscansharedtargets$TF)
    rownames(notsp) = unique(sp_interest_format_targetscansharedtargets$Target)
    
    for (i in 1:nrow(sp_interest_format_targetscansharedtargets))
      if (sp_interest_format_targetscansharedtargets[i, 'SignedQuantile'] > 0){
        notsp[sp_interest_format_targetscansharedtargets[i,'Target'], sp_interest_format_targetscansharedtargets[i,'TF']] = 1
      } else if (sp_interest_format_targetscansharedtargets[i, 'SignedQuantile'] < 0){
        notsp[sp_interest_format_targetscansharedtargets[i,'Target'], sp_interest_format_targetscansharedtargets[i,'TF']] = -1
      }
    
    notsp_new <- notsp[ order(row.names(notsp)), ]
    notsp_new <- notsp[ , order(colnames(notsp))]
    
    
    write.table(notsp_new, paste0(outputfolder, '/ms', ms, '/bias', bias, '/sharedbyMorethan', num, 'netsFrom5seeds_pos/drcombat_ts_ct02_bias', bias, tfaOpt, rm, '_cut01_sharedbyMorethan', num, 'nets_pos_targetscansharedtargets.tsv'), sep = '\t', quote = F, col.names=NA)
    
    
    notsp_new <- notsp[ order(row.names(notsp)), ]
    write.table(notsp_new, paste0(outputfolder, '/ms', ms, '/bias', bias, '/sharedbyMorethan', num, 'netsFrom5seeds_pos/drcombat_ts_ct02_bias', bias, tfaOpt, rm, '_cut01_sharedbyMorethan', num, 'nets_pos_targetscansharedtargets_test.tsv'), sep = '\t', quote = F, col.names=NA)
    
}
    
  ts_ct02_b = read.table('/workdir/hz543/projects/Inferelator/mirna_gene_naiveOnly_noSP8_wt_prior/inputs/prior/ts_ct02_contextscore.tsv')
  #num = 2
  
  for (num in 1:4){
    sp_interest_names = names(res_table_pos[res_table_pos > num])
    
    sp_interest = sp_all_pos[sp_all_pos$`TF->Target,direction` %in% sp_interest_names, ]
    
    sp_interest_less = sp_interest[, c('TF->Target,direction', 'TF.x')]
    sp_interest_less$TF = gsub('->.*', '', sp_interest_less$`TF->Target,direction`)
    sp_interest_less$Target = gsub(',.*', '', gsub('.*->', '', sp_interest_less$`TF->Target,direction`))
    sp_interest_less$direction = gsub('.*,', '', sp_interest_less$`TF->Target,direction`)
    sp_interest_less$SignedQuantile = NA
    for (i in 1:nrow(sp_interest_less)){
      if (sp_interest_less[i,'direction'] == 'pos'){
        sp_interest_less[i,'SignedQuantile'] = 1
      } else if (sp_interest_less[i,'direction'] == 'neg'){
        sp_interest_less[i,'SignedQuantile'] = -1
      }
    }
    sp_interest_format = sp_interest_less[, c('TF', 'Target', 'SignedQuantile')]
    print(num)
    print(nrow(sp_interest_format))
    
    print(table(unique(sp_interest_format$TF) %in% colnames(ts_ct02_b))) # all true - all TFs(miRNAs) are listed in targetscan
    print(table(unique(sp_interest_format$Target) %in% rownames(ts_ct02_b))) # not all true - not all targets are listed in targetscan
    
    sp_interest_format_targetscan = sp_interest_format[sp_interest_format$Target %in% rownames(ts_ct02_b),]
    print(table(unique(sp_interest_format_targetscan$Target) %in% rownames(ts_ct02_b))) # all true 
    
    shared_targets = c()
    for (i in 1:nrow(sp_interest_format_targetscan)){
      if (ts_ct02_b[sp_interest_format_targetscan[i,'Target'], sp_interest_format_targetscan[i,'TF']] != 0){
        shared_targets = c(shared_targets, i)
      }
    }
    sp_interest_format_targetscansharedtargets = sp_interest_format_targetscan[shared_targets, ]
    write.table(sp_interest_format_targetscansharedtargets, paste0(outputfolder, '/ms', ms, '/bias', bias, '/sharedbyMorethan', num, 'netsFrom5seeds_pos/drcombat_ts_ct02_bias', bias, tfaOpt, rm, '_cut01_sharedbyMorethan', num, 'nets_pos_targetscansharednopcttargets_sp.tsv'), sep = '\t', quote = F, row.names = F)  
    
    
    notsp = data.frame(matrix(0, ncol = length(unique(sp_interest_format_targetscansharedtargets$TF)), nrow = length(unique(sp_interest_format_targetscansharedtargets$Target))))
    colnames(notsp) = unique(sp_interest_format_targetscansharedtargets$TF)
    rownames(notsp) = unique(sp_interest_format_targetscansharedtargets$Target)
    
    for (i in 1:nrow(sp_interest_format_targetscansharedtargets))
      if (sp_interest_format_targetscansharedtargets[i, 'SignedQuantile'] > 0){
        notsp[sp_interest_format_targetscansharedtargets[i,'Target'], sp_interest_format_targetscansharedtargets[i,'TF']] = 1
      } else if (sp_interest_format_targetscansharedtargets[i, 'SignedQuantile'] < 0){
        notsp[sp_interest_format_targetscansharedtargets[i,'Target'], sp_interest_format_targetscansharedtargets[i,'TF']] = -1
      }
    
    notsp_new <- notsp[ order(row.names(notsp)), ]
    notsp_new <- notsp[ , order(colnames(notsp))]
    
    
    write.table(notsp_new, paste0(outputfolder, '/ms', ms, '/bias', bias, '/sharedbyMorethan', num, 'netsFrom5seeds_pos/drcombat_ts_ct02_bias', bias, tfaOpt, rm, '_cut01_sharedbyMorethan', num, 'nets_pos_targetscansharednopcttargets.tsv'), sep = '\t', quote = F, col.names=NA)
    
    
    notsp_new <- notsp[ order(row.names(notsp)), ]
    write.table(notsp_new, paste0(outputfolder, '/ms', ms, '/bias', bias, '/sharedbyMorethan', num, 'netsFrom5seeds_pos/drcombat_ts_ct02_bias', bias, tfaOpt, rm, '_cut01_sharedbyMorethan', num, 'nets_pos_targetscansharednopcttargets_test.tsv'), sep = '\t', quote = F, col.names=NA)
    
}

}


for (ms in c(2)) 
  for (rm in c('')) 
    for (bias in c(5))
      for (tfaOpt in c('_TFmRNA')){ 
        print(paste0('ms', ms, ', bias', bias, tfaOpt, rm))
        getsharednets(ms, bias, tfaOpt, rm)
      }

