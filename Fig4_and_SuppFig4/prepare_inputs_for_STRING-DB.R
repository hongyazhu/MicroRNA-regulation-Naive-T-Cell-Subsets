# get miR-Inf predicted miR-29 targets, conserved
net = read.table('analyze_network/ms2/bias5/prioritizeTargets/pos_regs_freqTFmRNAmorethan4.tsv', header = T)
net_mir29_con = net[net$miRNA_name == 'miR-29abc-3p' & net$ts %in% 'TRUE', c('miRNA_name', 'target')] # 
net_mir29_con$SignedQuantile = 1
write.table(net_mir29_con, 'make_figures/net_mir29_con.tsv', sep = '\t', quote = F, row.names = F)

net = read.table('analyze_network/ms2/bias5/prioritizeTargets/pos_regs_freqTFmRNAmorethan4.tsv', header = T)
net_let7_con = net[net$miRNA_name == 'let-7abcdefgik,98-5p,miR-1961' & net$ts %in% 'TRUE', c('miRNA_name', 'target')] #  & net$ts %in% 'TRUE'
net_let7_con$SignedQuantile = 1
write.table(net_let7_con, 'make_figures/ct02_ms2_4ormore/net_let7_con.tsv', sep = '\t', quote = F, row.names = F)


# background in STRING-DB: inputs/exonIntron_matrix/targGene_match_noNA_noInf_noall0_10exon_tsmirnasexpressed1000_ct02.txt 
# (generated in Fig3_and_SuppFig3/0_inputs_for_miR-Inf/generating_targGene_file.R)
