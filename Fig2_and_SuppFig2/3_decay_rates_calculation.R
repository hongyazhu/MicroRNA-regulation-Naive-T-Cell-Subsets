
library(ggplot2)
library(ggfortify)
suppressMessages(library(DESeq2))

# https://gist.github.com/slowkow/6e34ccb4d1311b8fe62e
tpm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(rate) * 1e6
}

# read in exon counts
exon_rc = read.table("carpRNAseq_counts_exon.txt", header = T, row.names = 1)
exon_feature_lengths = exon_rc$Length
exon_rc = exon_rc[,6:ncol(exon_rc)]
colnames(exon_rc) = c('Empty_rep1_batch1', 'Empty_rep2_batch1', 'miR.1_rep1_batch1', 'miR.1_rep2_batch1', 'miR.122_rep1_batch1', 'miR.122_rep2_batch1',
                      'Empty_rep1_batch2', 'Empty_rep2_batch2', 'Empty_rep3_batch2', 'miR.133a_rep1_batch2', 'miR.133a_rep2_batch2', 'miR.133a_rep3_batch2',
                      'miR.155_rep1_batch2', 'miR.155_rep2_batch2', 'miR.155_rep3_batch2', 'miR.302a_rep1_batch2', 'miR.302a_rep2_batch2', 'miR.302a_rep3_batch2',
                      'miR.372_rep1_batch2', 'miR.372_rep2_batch2', 'miR.372_rep3_batch2', 'miR.373_rep1_batch2', 'miR.373_rep2_batch2', 'miR.373_rep3_batch2')


# read in gene body counts
gb_rc = read.table("carpRNAseq_counts_gb.txt", header = T, row.names = 1)
gb_feature_lengths = gb_rc$Length
gb_rc = gb_rc[,6:ncol(gb_rc)]
colnames(gb_rc) = c('Empty_rep1_batch1', 'Empty_rep2_batch1', 'miR.1_rep1_batch1', 'miR.1_rep2_batch1', 'miR.122_rep1_batch1', 'miR.122_rep2_batch1',
                    'Empty_rep1_batch2', 'Empty_rep2_batch2', 'Empty_rep3_batch2', 'miR.133a_rep1_batch2', 'miR.133a_rep2_batch2', 'miR.133a_rep3_batch2',
                    'miR.155_rep1_batch2', 'miR.155_rep2_batch2', 'miR.155_rep3_batch2', 'miR.302a_rep1_batch2', 'miR.302a_rep2_batch2', 'miR.302a_rep3_batch2',
                    'miR.372_rep1_batch2', 'miR.372_rep2_batch2', 'miR.372_rep3_batch2', 'miR.373_rep1_batch2', 'miR.373_rep2_batch2', 'miR.373_rep3_batch2')

# calculating intron counts
intron_rc <- gb_rc - exon_rc
intron_feature_lengths = gb_feature_lengths - exon_feature_lengths

intron_rc_noNeg <- intron_rc[rowSums(intron_rc < 0) == 0,]
exon_rc_noNeg <- exon_rc[rowSums(intron_rc < 0) == 0,]

intron_length_noNeg = intron_feature_lengths[rowSums(intron_rc < 0) == 0]
exon_length_noNeg = exon_feature_lengths[rowSums(intron_rc < 0) == 0]

intron_rc_wIntron = intron_rc_noNeg[intron_length_noNeg != 0, ]
exon_rc_wIntron = exon_rc_noNeg[intron_length_noNeg != 0, ]

intron_length_wIntron = intron_length_noNeg[intron_length_noNeg != 0]
exon_length_wIntron = exon_length_noNeg[intron_length_noNeg != 0]

# normalize to tpm
tpm_intron_wIntron <- apply(intron_rc_wIntron, 2, function(x) tpm(x, intron_length_wIntron))
tpm_exon_wIntron <- apply(exon_rc_wIntron, 2, function(x) tpm(x, exon_length_wIntron))

# calculating decay rates with exon and intron tpm
decay_rates = data.frame(matrix(nrow = nrow(tpm_exon_wIntron), ncol = ncol(tpm_exon_wIntron)))
colnames(decay_rates) = colnames(tpm_exon_wIntron)
rownames(decay_rates) = rownames(tpm_exon_wIntron)
for (i in 1:nrow(decay_rates)){
    for (j in 1:ncol(decay_rates)){
        decay_rates[i,j] = tpm_intron_wIntron[i,j] / tpm_exon_wIntron[i,j]
    }
}
write.csv(decay_rates, "decay_rates.csv", quote = F)
                          
