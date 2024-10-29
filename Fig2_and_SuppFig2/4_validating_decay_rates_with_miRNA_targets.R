# Evaluate decay rates by comparing decay rates of miRNA targets and background in control and miRNA-overexpressed samples.
# Data from Patel et al., 2020
# Fig 2C-D, Supp Fig 2D-F

library(VennDiagram)
library(edgeR)
library(reshape2)
library(ggplot2)

### Prep: find targets from Targetscan 

# from Targetscan, find miRNAs' seed sequence (used to find targets)
mir_family_info <- read.table("/path/to/targetscan/human/miR_Family_Info.txt", head = T, fill = T, sep = "\t")
hsa_mir_family_info <- mir_family_info[substr(mir_family_info[,'MiRBase.ID'] ,1 ,3) == "hsa", ]
hsa_con_mir_family_info <- hsa_mir_family_info[hsa_mir_family_info[,'Family.Conservation.'] >= 2, ]
seed_con <- hsa_con_mir_family_info[,c('Seed.m8', 'MiRBase.ID')]
mirna_seeds_con <- aggregate(seed_con[,2], by = list(seed = seed_con$Seed), FUN = paste)
colnames(mirna_seeds_con) <- c("seed", "miRNA")

# find targets (with transcript ID) (filtered by total context score < -0.2)
targets_human <- read.table("/path/to/targetscan/human/Targetscan_human_filtered_genes.txt")
colnames(targets_human) <- c("transcript_id_ver", "gene_name","seed", "mirna_name")
targets_human <- targets_human[,c("gene_name", "seed")]

# find weak targets (to be excluded in the analysis) (filted by total context score > -0.2)
weak_targets <- read.table("/path/to/data/targetscan/human/Targetscan_human_weak_targets_genesNames.txt")
colnames(weak_targets) <- c("gene_name", "seed")

# file for converting ensembl id and gene name
id2name = read.table("/path/to/human/gencode_annotation/ensemble_ID2gene_name.txt")
colnames(id2name) = c('id', 'gene_name')


### Read decay rate data
decay_rates = read.table("decay_rates.csv", sep = ",", header = T, row.names = 1)
rownames(decay_rates) = substr(rownames(decay_rates), 1, 15)



# Function: plot figures for decay rates of miRNA targets predicted by CARP, RNAseq and Targetscan, log2 decay rates for boxplots
### input
# miRNA, e.g. 'm1'
# miRNA.name_, e.g. 'miRNA.name_'
# seed, e.g. 'GGAAUGU', 'UUGGUCC|UGGUCCC'
# batch, e.g. 'batch1'
### output
# venn diagram for targets shared by the three methods
# boxplots comparing decay rates for predicted miRNA targets (and background set of genes) in relevant sample and ctrl, log2 decay rates
# lineplots comparing logFC of decay rates for predicted targets and background genes between relevant sample and ctrl

decay_rates_figures <- function(miRNA, miRNA.name_, seed, batch){

    #############################################################
    ################### find relevant samples ###################
    #############################################################

    sample_name_batch = grep(batch, colnames(decay_rates), value=TRUE)
    sample_name_rel = c(grep('Empty_', sample_name_batch, value=TRUE), grep(miRNA.name_, sample_name_batch, value=TRUE))


    #############################################################
    ############ find targets identified by methods #############
    #############################################################

    # targets - criteria from Ravi

    # genes that are not lowly expressed - only those are included in the analysis
    #  exon cpm, relevant, not low
    exon_rc = read.table("carpRNAseq_counts_exon.txt", header = T, row.names = 1)
    exon_feature_lengths = exon_rc$Length
    exon_rc = exon_rc[,6:ncol(exon_rc)]
    colnames(exon_rc) = c('Empty_rep1_batch1', 'Empty_rep2_batch1', 'miR.1_rep1_batch1', 'miR.1_rep2_batch1', 'miR.122_rep1_batch1', 'miR.122_rep2_batch1',
                          'Empty_rep1_batch2', 'Empty_rep2_batch2', 'Empty_rep3_batch2', 'miR.133a_rep1_batch2', 'miR.133a_rep2_batch2', 'miR.133a_rep3_batch2',
                          'miR.155_rep1_batch2', 'miR.155_rep2_batch2', 'miR.155_rep3_batch2', 'miR.302a_rep1_batch2', 'miR.302a_rep2_batch2', 'miR.302a_rep3_batch2',
                          'miR.372_rep1_batch2', 'miR.372_rep2_batch2', 'miR.372_rep3_batch2', 'miR.373_rep1_batch2', 'miR.373_rep2_batch2', 'miR.373_rep3_batch2')
    exon_rc_rel = exon_rc[,sample_name_rel]
    exon_cpm = cpm(exon_rc_rel)
    exon_cpm_noLow = exon_cpm[rowSums(exon_cpm < 1) == 0 , , drop = FALSE] # 12874 genes
    exon_cpm_noLow = rownames(exon_cpm_noLow)

    # gene body cpm, relevant, not low
    gb_rc = read.table("carpRNAseq_counts_gb.txt", header = T, row.names = 1)
    gb_feature_lengths = gb_rc$Length
    gb_rc = gb_rc[,6:ncol(gb_rc)]
    colnames(gb_rc) = c('Empty_rep1_batch1', 'Empty_rep2_batch1', 'miR.1_rep1_batch1', 'miR.1_rep2_batch1', 'miR.122_rep1_batch1', 'miR.122_rep2_batch1',
                        'Empty_rep1_batch2', 'Empty_rep2_batch2', 'Empty_rep3_batch2', 'miR.133a_rep1_batch2', 'miR.133a_rep2_batch2', 'miR.133a_rep3_batch2',
                        'miR.155_rep1_batch2', 'miR.155_rep2_batch2', 'miR.155_rep3_batch2', 'miR.302a_rep1_batch2', 'miR.302a_rep2_batch2', 'miR.302a_rep3_batch2',
                        'miR.372_rep1_batch2', 'miR.372_rep2_batch2', 'miR.372_rep3_batch2', 'miR.373_rep1_batch2', 'miR.373_rep2_batch2', 'miR.373_rep3_batch2')
    gb_rc_rel = gb_rc[,sample_name_rel]
    gb_cpm = cpm(gb_rc_rel)
    gb_cpm_noLow = gb_cpm[rowSums(gb_cpm < 1) == 0 , , drop = FALSE] # 12874 genes
    gb_cpm_noLow = rownames(gb_cpm_noLow)


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


    ### CARP targets (file downloaded from GSE142895)
    # CARP edgeR: q-value < 0.05 and log2(miR/Control) < –0.2
    miRNA_CARP_edgeR = read.table(paste0("GSE142895_edgeR.CARP.", miRNA, ".vs.Empty.tsv"), header = T)
    miRNA_CARP = miRNA_CARP_edgeR[miRNA_CARP_edgeR$qValue < 0.05 & miRNA_CARP_edgeR$logFC < -0.2, ]
    miRNA_CARP = substr(miRNA_CARP$Geneid, 1, 15)

    miRNA_CARP_name = id2name[id2name$id %in% miRNA_CARP, 'gene_name'] # 288 -> 285

    # include only genes that are not lowly expressed
    miRNA_CARP_notLow = intersect(miRNA_CARP_name, exon_cpm_noLow)
    miRNA_CARP_notLow = intersect(miRNA_CARP_notLow, gb_cpm_noLow)

    # remove weak targets
    miRNA_CARP_noWeak = setdiff(miRNA_CARP_notLow, miRNA_weak_ts)


    ### RNAseq targets (file downloaded from GSE142895)
    # RNAseq edgeR: q-value < 0.05 and log2(miR/Control) < –0.2
    miRNA_RNA_edgeR = read.table(paste0("GSE142895_edgeR.rnaSeq.", miRNA, ".vs.Empty.tsv"), header = T)
    miRNA_RNA = miRNA_RNA_edgeR[miRNA_RNA_edgeR$qValue < 0.05 & miRNA_RNA_edgeR$logFC < -0.2, ]
    miRNA_RNA = substr(miRNA_RNA$Geneid, 1, 15)

    miRNA_RNA_name = id2name[id2name$id %in% miRNA_RNA, 'gene_name'] # 1748 -> 1714

    # include only genes that are not lowly expressed
    miRNA_RNA_notLow = intersect(miRNA_RNA_name, exon_cpm_noLow)
    miRNA_RNA_notLow = intersect(miRNA_RNA_notLow, gb_cpm_noLow)

    # remove weak targets
    miRNA_RNA_noWeak = setdiff(miRNA_RNA_notLow, miRNA_weak_ts)

    ### Targetscan targets
    # context++ score < –0.2
    if (!grepl('\\|', seed)){
        miRNA_ts <- targets_human[targets_human$seed == seed,]$gene_name
    } else {
        seed1 = strsplit(seed,split='|', fixed=T)[[1]][1]
        seed2 = strsplit(seed,split='|', fixed=T)[[1]][2]
        miRNA_ts1 <- targets_human[targets_human$seed == seed1,]$gene_name
        miRNA_ts2 <- targets_human[targets_human$seed == seed2,]$gene_name
        miRNA_ts = union(miRNA_ts1, miRNA_ts2)
    }


    #miRNA_ts <- targets_human_geneID[targets_human_geneID$seed == seed,]$gene_id

    # include only genes that are not lowly expressed
    miRNA_ts_notLow = intersect(miRNA_ts, exon_cpm_noLow)
    miRNA_ts_notLow = intersect(miRNA_ts_notLow, gb_cpm_noLow)

    # remove weak targets
    # miRNA_ts_noWeak = setdiff(miRNA_ts_notLow, miRNA_weak_ts)


    #############################################################
    ################# find decay rates of targets ###############
    #############################################################

    # decay rates relevant to this miRNA
    decay_rates_rel = decay_rates[, sample_name_rel]

    # remove genes with NAs and Inf decay rates (dr)
    decay_rates_rel_noNA <- na.omit(decay_rates_rel)
    decay_rates_rel_noNA_noInf <- decay_rates_rel_noNA[!is.infinite(rowSums(decay_rates_rel_noNA)),]

    ### select decay rates with targets

    miRNA_CARP_dr = decay_rates_rel_noNA_noInf[intersect(miRNA_CARP_noWeak, rownames(decay_rates_rel_noNA_noInf)), ]

    miRNA_RNA_dr = decay_rates_rel_noNA_noInf[intersect(miRNA_RNA_noWeak, rownames(decay_rates_rel_noNA_noInf)), ]

    miRNA_ts_dr = decay_rates_rel_noNA_noInf[intersect(miRNA_ts_notLow, rownames(decay_rates_rel_noNA_noInf)), ]


    ### select decay rates with background genes
    # include only genes that are not lowly expressed
    bg_notLow = intersect(rownames(decay_rates_rel_noNA_noInf), exon_cpm_noLow)
    bg_notLow = intersect(bg_notLow, gb_cpm_noLow)

    # remove weak targets
    bg_noWeak = setdiff(bg_notLow, miRNA_weak_ts)

    # remove all targets
    bg_noTargets = setdiff(bg_noWeak, miRNA_CARP_dr)
    bg_noTargets = setdiff(bg_noTargets, miRNA_RNA_dr)
    bg_noTargets = setdiff(bg_noTargets, miRNA_ts_dr)
    miRNA_bg_dr = decay_rates_rel_noNA_noInf[intersect(bg_noTargets, rownames(decay_rates_rel_noNA_noInf)), ]

    #############################################################
    ######################## plot boxplots ######################
    #############################################################

    miRNA_CARP_dr_toMelt = log2(miRNA_CARP_dr)
    miRNA_RNA_dr_toMelt = log2(miRNA_RNA_dr)
    miRNA_ts_dr_toMelt = log2(miRNA_ts_dr)
    miRNA_bg_dr_toMelt = log2(miRNA_bg_dr)

    miRNA_CARP_dr_toMelt$method = 'CARP'
    miRNA_RNA_dr_toMelt$method = 'RNAseq'
    miRNA_ts_dr_toMelt$method = 'Targetscan'
    miRNA_bg_dr_toMelt$method = 'Background'

    miRNA_CARP_dr_melt = melt(miRNA_CARP_dr_toMelt)
    miRNA_RNA_dr_melt = melt(miRNA_RNA_dr_toMelt)
    miRNA_ts_dr_melt = melt(miRNA_ts_dr_toMelt)
    miRNA_bg_dr_melt = melt(miRNA_bg_dr_toMelt)

    miRNA_dr_melt = rbind(rbind(miRNA_CARP_dr_melt, miRNA_RNA_dr_melt), rbind(miRNA_ts_dr_melt, miRNA_bg_dr_melt))
    miRNA_dr_melt_finite <- miRNA_dr_melt[!is.infinite(miRNA_dr_melt$value),]
    print(nrow(miRNA_dr_melt))
    print(nrow(miRNA_dr_melt_finite))
    
    # find number of targets
    CARP_ntargets = nrow(miRNA_CARP_dr)
    RNA_ntargets = nrow(miRNA_RNA_dr)
    ts_ntargets = nrow(miRNA_ts_dr)
    bg_ntargets = nrow(miRNA_bg_dr)

    # Change box plot colors by groups
    if (batch == "batch1"){
        options(repr.plot.width = 7, repr.plot.height = 5)
        boxp <- ggplot(miRNA_dr_melt_finite, aes(x=method, y=value, fill=variable)) +
            geom_boxplot() +
            xlab("Method") +
            ylab("Log2 decay rates") +
            theme_bw() +
            scale_fill_manual(values=c("#999999", "#999999", "#E69F00", "#E69F00")) +
            theme(text=element_text(size=15),
                plot.title = element_text(hjust = 0.5)) +
            ggtitle(paste0(miRNA.name_, "Targets Decay Rates")) +
            scale_x_discrete(limits = c('Background', 'CARP', 'RNAseq', 'Targetscan'),
                             labels = paste0(c('Background', 'CARP', 'RNAseq', 'Targetscan'), "\nn=", c(bg_ntargets, CARP_ntargets, RNA_ntargets, ts_ntargets)))
    } else {
        options(repr.plot.width = 7, repr.plot.height = 5)
        boxp <- ggplot(miRNA_dr_melt_finite, aes(x=method, y=value, fill=variable)) +
            geom_boxplot() +
            xlab("Method") +
            ylab("Log2 decay rates") +
            theme_bw() +
            scale_fill_manual(values=c("#999999", "#999999", "#999999", "#E69F00", "#E69F00", "#E69F00")) +
            theme(text=element_text(size=15),
                plot.title = element_text(hjust = 0.5)) +
            ggtitle(paste0(miRNA.name_, "Targets Decay Rates")) +
            scale_x_discrete(limits = c('Background', 'CARP', 'RNAseq', 'Targetscan'),
                             labels = paste0(c('Background', 'CARP', 'RNAseq', 'Targetscan'), "\nn=", c(bg_ntargets, CARP_ntargets, RNA_ntargets, ts_ntargets)))
        
    }

    ggsave(
        paste0("compare_tp_methods/boxplots_log2/", miRNA.name_, "targets_dr.pdf"),
        plot = boxp,
        device = "pdf",
        width = 7,
        height = 5,
        dpi = 300
        )

    #############################################################
    ########################### plot fc #########################
    #############################################################

    miRNA_CARP_dr_fc = miRNA_CARP_dr
    miRNA_RNA_dr_fc = miRNA_RNA_dr
    miRNA_ts_dr_fc = miRNA_ts_dr
    miRNA_bg_dr_fc = miRNA_bg_dr

    ctrl_sample_rel = grep("Empty", sample_name_rel, value=TRUE)
    miRNA_sample_rel = sample_name_rel[!grepl("Empty", sample_name_rel)]

    miRNA_CARP_dr_fc$fc = log2(rowMeans(miRNA_CARP_dr_fc[,miRNA_sample_rel])/rowMeans(miRNA_CARP_dr_fc[,ctrl_sample_rel]))
    miRNA_RNA_dr_fc$fc = log2(rowMeans(miRNA_RNA_dr_fc[,miRNA_sample_rel])/rowMeans(miRNA_RNA_dr_fc[,ctrl_sample_rel]))
    miRNA_ts_dr_fc$fc = log2(rowMeans(miRNA_ts_dr_fc[,miRNA_sample_rel])/rowMeans(miRNA_ts_dr_fc[,ctrl_sample_rel]))
    miRNA_bg_dr_fc$fc = log2(rowMeans(miRNA_bg_dr_fc[,miRNA_sample_rel])/rowMeans(miRNA_bg_dr_fc[,ctrl_sample_rel]))
        
    wilcox_test_CARP <- wilcox.test(miRNA_CARP_dr_fc$fc, miRNA_bg_dr_fc$fc)
    wilcox_test_RNA <- wilcox.test(miRNA_RNA_dr_fc$fc, miRNA_bg_dr_fc$fc)
    wilcox_test_ts <- wilcox.test(miRNA_ts_dr_fc$fc, miRNA_bg_dr_fc$fc)
    wilcox_test_bg <- wilcox.test(miRNA_bg_dr_fc$fc, miRNA_bg_dr_fc$fc)

    options(repr.plot.width = 6, repr.plot.height = 6)
    pdf(paste0("compare_tp_methods/lineplots/", miRNA.name_, "targets_fc.pdf"), width = 8, height = 8) 
    
    plot(as.list(environment(ecdf(miRNA_CARP_dr_fc$fc))), pch = ".", col = "white", xlim=c(-1,1), 
        xlab = paste0("Log Fold Change: ", miRNA.name_, "vs_", "Ctrl"), ylab = "Density", 
        cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
    lines(as.list(environment(ecdf(miRNA_CARP_dr_fc$fc))), col = '#cf1578', pch = ".", lwd = 3)
    lines(as.list(environment(ecdf(miRNA_RNA_dr_fc$fc))), col = '#e8d21d', pch = ".", lwd = 3)
    lines(as.list(environment(ecdf(miRNA_ts_dr_fc$fc))), col = '#039fbe', pch = ".", lwd = 3)
    lines(as.list(environment(ecdf(miRNA_bg_dr_fc$fc))), col = 'black', lwd = 3)
    abline(v=0, h=0.5, col = c("gray","gray"), lty = c(3,3))  
    legend(-1.05, 1, c(paste0("Background (n=", nrow(miRNA_bg_dr_fc), ")"),
                    paste0("CARP (n=", nrow(miRNA_CARP_dr_fc), ", p=", formatC(wilcox_test_CARP$p.value, format = "e", digits = 2), ")"), 
                    paste0("RNAseq (n=", nrow(miRNA_RNA_dr_fc), ", p=", formatC(wilcox_test_RNA$p.value, format = "e", digits = 2), ")"), 
                    paste0("Targetscan (n=", nrow(miRNA_ts_dr_fc), ", p=", formatC(wilcox_test_ts$p.value, format = "e", digits = 2), ")")), 
                col = c("black", "#cf1578", "#e8d21d", "#039fbe"), 
                lty=c(1,1), 
                cex = 0.9, 
                lwd = 1.5,
                pt.cex = 10)
    title(main = paste0(miRNA.name_, "Targets Decay Rates Fold Change"), cex.main = 1.5)

    dev.off()
}

decay_rates_figures(miRNA = "m1", miRNA.name_ = "miR.1_", seed = 'GGAAUGU', batch = 'batch1')
decay_rates_figures(miRNA = "m122", miRNA.name_ = "miR.122_", seed = 'GGAGUGU', batch = 'batch1')

decay_rates_figures(miRNA = "m133a", miRNA.name_ = "miR.133a_", seed = 'UUGGUCC|UGGUCCC', batch = 'batch2')
decay_rates_figures(miRNA = "m155", miRNA.name_ = "miR.155_", seed = 'UAAUGCU', batch = 'batch2')
decay_rates_figures(miRNA = "m302a", miRNA.name_ = "miR.302a_", seed = 'AAGUGCU', batch = 'batch2')
decay_rates_figures(miRNA = "m372", miRNA.name_ = "miR.372_", seed = 'AAGUGCU', batch = 'batch2')
decay_rates_figures(miRNA = "m373", miRNA.name_ = "miR.373_", seed = 'AAGUGCU', batch = 'batch2')
