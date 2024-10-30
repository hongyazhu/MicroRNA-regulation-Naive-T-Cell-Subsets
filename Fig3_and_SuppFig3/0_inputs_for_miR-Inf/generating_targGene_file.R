# generating target gene file, including genes that have decay rates estimated, and bear miRNA target sites based on TargetScan.

targGeneAll = read.table(paste0('inputs/exonIntron_matrix/targGene_match_noNA_noInf_noall0_10exon.txt')) 
# generated in Fig2_and_SuppFig2/5_decay_rates_for_naive_T_cell_datasets.R

tsrel_ct02 = read.table('inputs/prior/tsmirnasexpressed1000_ct02_contextscore0to1.tsv')
# generated in generating_prior_from_TargetScan.R

targGene_tsrel = intersect(rownames(tsrel_ct02), targGeneAll$V1)
write.table(targGene_tsrel, paste0("exonIntron_matrix/targGene_match_noNA_noInf_noall0_10exon_tsmirnasexpressed1000_ct02.txt"), row.names = FALSE, col.names=FALSE, sep = "\t", quote = F)
