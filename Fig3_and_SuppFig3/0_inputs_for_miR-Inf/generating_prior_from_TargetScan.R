
outputfolder = 'inputs/prior/'

targetscan_reg = read.table("/path/to/targetscan/mouse/Targetscan_mouse_ct02.txt") # TargetScan predictions, filtered by total context score <- 0.2
colnames(targetscan_reg) = c('transcript_id', 'gene_name', 'seed', 'mirna', 'total_context_score', 'cumulative_weighted_context_score', 'aggregate_pct')
targetscan_reg = targetscan_reg[,c('gene_name', 'seed', 'total_context_score')]

mir_family_info <- read.table("/path/to/targetscan/mouse/miR_Family_Info.txt", head = T, fill = T, sep = "\t")
mmu_mir_family_info <- mir_family_info[substr(mir_family_info[,'MiRBase.ID'] ,1 ,3) == "mmu", ]
mmu_con_mir_family_info <- mmu_mir_family_info[mmu_mir_family_info[,'Family.Conservation.'] >= 2, ]
seed_con <- mmu_con_mir_family_info[,c('Seed.m8', 'MiRBase.ID')]
mirna_seeds_con <- aggregate(seed_con[,2], by = list(seed = seed_con$Seed), FUN = paste)
colnames(mirna_seeds_con) <- c("seed", "miRNA")

targetscan_reg_seed = merge(targetscan_reg, mirna_seeds_con, by='seed')
targetscan_reg_seed$miRNA <- vapply(targetscan_reg_seed$miRNA, paste, collapse = ";", character(1L))

prior_mat = data.frame(matrix(0, nrow = length(unique(targetscan_reg_seed$gene_name)), ncol = length(unique(targetscan_reg_seed$seed))))

rownames(prior_mat) = unique(targetscan_reg_seed$gene_name)
colnames(prior_mat) = unique(targetscan_reg_seed$seed)

for (i in 1:nrow(targetscan_reg_seed)){
  mirna_tmp = targetscan_reg_seed[i, 'seed']
  target_tmp = targetscan_reg_seed[i, 'gene_name']
  prior_mat[target_tmp, mirna_tmp] = targetscan_reg_seed[i, 'total_context_score']
}

write.table(prior_mat, paste0(outputfolder, "ts_ct02_contextscore.tsv"), quote = F, sep = "\t", col.names=NA)



# generate prior from targetscan retaining only regulators that are expressed
ts = read.table('inputs/prior/ts_ct02_contextscore.tsv')
expressed_minras = read.table('inputs/mirna_DE/mirnas_expressed1000.txt')$V1
ts_rel = ts[, colnames(ts) %in% expressed_minras]
ts_rel = ts_rel[rowSums(ts_rel[])!=0,]
ts_rel = -ts_rel
ts_rel = ts_rel/max(ts_rel)
write.table(ts_rel, "/workdir/hz543/projects/Inferelator/mirna_gene_naiveOnly_noSP8_wt_prior/inputs/prior/tsmirnasexpressed1000_ct02_contextscore0to1.tsv", quote = F, sep = "\t", col.names=NA)
# run generating_targGene_file.R to get targGene file (targGene_match_noNA_noInf_noall0_10exon_tsmirnasexpressed1000_ct02.txt)



# scale values in prior matrix (which were total context scores) to 0 to 1
ts = read.table('inputs/prior/ts_ct02_contextscore.tsv')
expressed_minras = read.table('inputs/mirna_DE/mirnas_expressed1000.txt')$V1
ts_rel = ts[, colnames(ts) %in% expressed_minras]
ts_rel = ts_rel[rowSums(ts_rel[])!=0,]
targGene = read.table('inputs/exonIntron_matrix/targGene_match_noNA_noInf_noall0_10exon_tsmirnasexpressed1000_ct02.txt')$V1
ts_rel = ts_rel[targGene, ]
ts_rel = -ts_rel
for (i in 1:nrow(ts_rel))
  for (j in 1:ncol(ts_rel))
    if (ts_rel[i,j] != 0)
      ts_rel[i,j] = (ts_rel[i,j] - min(ts_rel[ts_rel > 0]))/max(ts_rel)
hist(as.matrix(ts_rel), xlab = 'log10(-context scores)', main = 'Histogram of context scores\n(miRNA TPM>1000 in >=1 samples\ncontext scores <-0.1)')
write.table(ts_rel, "inputs/prior/tsmirnasexpressed1000_ct02_contextscore0to1_targGeneonly.tsv", quote = F, sep = "\t", col.names=NA)

