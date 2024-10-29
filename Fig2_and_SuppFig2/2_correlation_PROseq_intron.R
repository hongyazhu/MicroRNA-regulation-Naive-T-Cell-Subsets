# Calculating correlation between PRO-seq and intron read counts
# Fig 2A, Supp Fig 2A-C

tpm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(rate) * 1e6
}


# read in PROseq data (Patel et al., 2020), downloaded from GEO
batch1_PROgb_rc = read.csv("GSE140365_batch1.genebody.read.count.tsv", header = T, sep = "\t") # 22020 genes
batch2_PROgb_rc = read.csv("GSE140365_batch2.genebody.read.count.tsv", header = T, sep = "\t") # 22033 genes

# record feature lengths
batch1_PRO_feature_lengths = batch1_PROgb_rc$Feature_length
batch2_PRO_feature_lengths = batch2_PROgb_rc$Feature_length

# drop the feature_length column
batch1_PROgb_rc = batch1_PROgb_rc[ , !(names(batch1_PROgb_rc) %in% c('Feature_length'))]
batch2_PROgb_rc = batch2_PROgb_rc[ , !(names(batch2_PROgb_rc) %in% c('Feature_length'))]
### feature lengths different for the two batches
# table(merge_PROgb_rc$Feature_length.x == merge_PROgb_rc$Feature_length.y)
#FALSE  TRUE 
# 8857 13121 

# take the gene ID column as row names
batch1_PROgb_rc <- data.frame(batch1_PROgb_rc, row.names = 1)
batch2_PROgb_rc <- data.frame(batch2_PROgb_rc, row.names = 1)

# normalization to tpm (previous normalization is wrong)
batch1_PROgb_rc_tpm <- apply(batch1_PROgb_rc, 2, function(x) tpm(x, batch1_PRO_feature_lengths))
batch2_PROgb_rc_tpm <- apply(batch2_PROgb_rc, 2, function(x) tpm(x, batch2_PRO_feature_lengths))

# merge two batches
merge_PROgb_rc_tpm = merge(batch1_PROgb_rc_tpm, batch2_PROgb_rc_tpm, by = 0) # 21978 genes

# take the gene ID column as row names
merge_PROgb_rc_tpm <- data.frame(merge_PROgb_rc_tpm, row.names = 1)
                        
                             
                             
                             
exon_rc = read.table("carpRNAseq_counts_exon.txt", header = T, row.names = 1)
exon_feature_lengths = exon_rc$Length
exon_rc = exon_rc[,6:ncol(exon_rc)]
colnames(exon_rc) = c('Empty_rep1_batch1', 'Empty_rep2_batch1', 'miR.1_rep1_batch1', 'miR.1_rep2_batch1', 'miR.122_rep1_batch1', 'miR.122_rep2_batch1',
                      'Empty_rep1_batch2', 'Empty_rep2_batch2', 'Empty_rep3_batch2', 'miR.133a_rep1_batch2', 'miR.133a_rep2_batch2', 'miR.133a_rep3_batch2',
                      'miR.155_rep1_batch2', 'miR.155_rep2_batch2', 'miR.155_rep3_batch2', 'miR.302a_rep1_batch2', 'miR.302a_rep2_batch2', 'miR.302a_rep3_batch2',
                      'miR.372_rep1_batch2', 'miR.372_rep2_batch2', 'miR.372_rep3_batch2', 'miR.373_rep1_batch2', 'miR.373_rep2_batch2', 'miR.373_rep3_batch2')



gb_rc = read.table("carpRNAseq_counts_gb.txt", header = T, row.names = 1)
gb_feature_lengths = gb_rc$Length
gb_rc = gb_rc[,6:ncol(gb_rc)]
colnames(gb_rc) = c('Empty_rep1_batch1', 'Empty_rep2_batch1', 'miR.1_rep1_batch1', 'miR.1_rep2_batch1', 'miR.122_rep1_batch1', 'miR.122_rep2_batch1',
                    'Empty_rep1_batch2', 'Empty_rep2_batch2', 'Empty_rep3_batch2', 'miR.133a_rep1_batch2', 'miR.133a_rep2_batch2', 'miR.133a_rep3_batch2',
                    'miR.155_rep1_batch2', 'miR.155_rep2_batch2', 'miR.155_rep3_batch2', 'miR.302a_rep1_batch2', 'miR.302a_rep2_batch2', 'miR.302a_rep3_batch2',
                    'miR.372_rep1_batch2', 'miR.372_rep2_batch2', 'miR.372_rep3_batch2', 'miR.373_rep1_batch2', 'miR.373_rep2_batch2', 'miR.373_rep3_batch2')

# getting intron read counts
intron_rc <- gb_rc - exon_rc
intron_feature_lengths = gb_feature_lengths - exon_feature_lengths

intron_rc_noNeg <- intron_rc[rowSums(intron_rc < 0) == 0,]
intron_length_noNeg = intron_feature_lengths[rowSums(intron_rc < 0) == 0]
intron_rc_wIntron = intron_rc_noNeg[intron_length_noNeg != 0, ]
intron_length_wIntron = intron_length_noNeg[intron_length_noNeg != 0]
tpm_intron_wIntron <- apply(intron_rc_wIntron, 2, function(x) tpm(x, intron_length_wIntron))


                                       
                            
tpm_intron_wIntron = as.data.frame(tpm_intron_wIntron)
tpm_intron_wIntron$gene_name = rownames(tpm_intron_wIntron)

rownames(merge_PROgb_rc_tpm) = substr(rownames(merge_PROgb_rc_tpm), 1, 15)

# file for converting ensembl id and gene name
id2name = read.table("/path/to/gencode_annotation/ensemble_ID2gene_name.txt")
colnames(id2name) = c('id', 'gene_name')

tpm_intron_wIntron_id = merge(tpm_intron_wIntron, id2name, by = 'gene_name')
rownames(tpm_intron_wIntron_id) = tpm_intron_wIntron_id$id
tpm_intron_wIntron_id$gene_name = NULL
tpm_intron_wIntron_id$id = NULL

intronPro = merge(tpm_intron_wIntron_id, merge_PROgb_rc_tpm, by = 0, suffixes = c(".intron",".pro"))
rownames(intronPro) = intronPro$Row.names
intronPro$Row.names = NULL


library(ggplot2)


sample_names = c('Empty_rep1_batch1', 'Empty_rep2_batch1', 'miR.1_rep1_batch1', 'miR.1_rep2_batch1', 'miR.122_rep1_batch1', 'miR.122_rep2_batch1',
                 'Empty_rep1_batch2', 'Empty_rep2_batch2', 'Empty_rep3_batch2', 'miR.133a_rep1_batch2', 'miR.133a_rep2_batch2', 'miR.133a_rep3_batch2',
                 'miR.155_rep1_batch2', 'miR.155_rep2_batch2', 'miR.155_rep3_batch2', 'miR.302a_rep1_batch2', 'miR.302a_rep2_batch2', 'miR.302a_rep3_batch2',
                 'miR.372_rep1_batch2', 'miR.372_rep2_batch2', 'miR.372_rep3_batch2', 'miR.373_rep1_batch2', 'miR.373_rep2_batch2', 'miR.373_rep3_batch2')

for(i in 1:length(sample_names)){
    sample_name = sample_names[i]
    
    intronPro_sample = intronPro[, paste0(sample_name, c('.intron', '.pro'))]

    intronPro_sample_no0 = intronPro_sample[apply(intronPro_sample, 1, function(row) all(row > 0 )),]
    intronPro_sample_no0_log2 = log2(intronPro_sample_no0)
    intronPro_sample_no0_log2 <- na.omit(intronPro_sample_no0_log2)
    intronPro_sample_no0_log2 <- intronPro_sample_no0_log2[!is.infinite(rowSums(intronPro_sample_no0_log2)),]



    x = intronPro_sample_no0_log2[,1]
    y = intronPro_sample_no0_log2[,2]
    m = lm(y ~ x)
    a = format(unname(coef(m)[1]), digits = 2)
    b = format(unname(coef(m)[2]), digits = 2)
    r2 = format(summary(m)$r.squared, digits = 3)

    p <- ggplot(intronPro_sample_no0_log2, aes(x = x, y = y)) + 
        geom_point(alpha = 0.1, color = 'green4') +
        xlab("log2 (Introns TPM)") +
        ylab("log2 (PROseq TPM)") +
        ggtitle(paste0(sample_name))  +
        geom_smooth(method = "lm", se=FALSE, color='black', formula = y ~ x, linetype = "longdash", size = 2) +
        theme_classic(base_size=16, base_family='ArialMT') +
        theme(plot.title = element_text(hjust = 0.5)) + 
        annotate("text", x = 0, y = 6, label = paste0("n = ", nrow(intronPro_sample_no0_log2)), size = 6) +
        annotate("text", x = 0, y = 8, label = paste0("y = ", a, " + ",b, "x , r^2 = ", r2), size = 6) 
                                                                         
    ggsave(
    paste0("correlation_intronPRO/", sample_name, ".pdf"),
    plot = p,
    device = "pdf",
    width = 6,
    height = 6,
    dpi = 300
    )


}


