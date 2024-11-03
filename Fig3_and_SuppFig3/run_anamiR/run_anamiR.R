library(anamiR)

## reading input
mirna = read.table('inputs/mirna_exp_matrix/counts_vst_blindT_combat_expressed1000.txt',
                   header = T)

mrna = read.table("make_figures/benchmark/anamiR/exon_countscombat_match.txt",
                  header = T)

pheno = data.frame(colnames(mirna), 
                                   c(rep('Wang 2016 (unpublished)', 7), rep('Wissink 2015', 8-2), rep('Clean CD8', 12), rep('Mir29Cre', 12-4), rep('CD8 Dev', 16-8)),# rep('Wang 2016', 2), ep('Clean Thy', 7), 
                                   c(#'Adult', 'Adult', 
                                     rep('Adult', 4), rep('Neo', 3),
                                     rep('Adult', 2), rep('Neo', 2), rep(c('Adult', 'Neo'), 2-1),
                                     rep('Adult', 12), 
                                     #rep('Adult', 7), 
                                     rep('Adult', 12-4),
                                     rep(c('Adult', 'Neo', 'Adult', 'Neo'), 2)),
                                   c(#'All', 'All', 
                                     'MP', 'MP', 'TN', 'TN', 'MP', 'MP', 'TN', 
                                     rep('All', 8-2),
                                     'TN', 'TN', 'TN', 'MP', 'MP', 'MP', 'TN', 'TN', 'TN', 'MP', 'MP', 'MP',
                                     #rep('All', 7),
                                     'MP', 'MP', 'TN', 'TN', 'MP', 'MP', 'TN', 'TN', # 'MP', 'MP', 'TN', 'TN', 
                                     rep('All', 16-8)),
                                   c(#rep('Spleen', 2), 
                                     rep('Spleen', 7), rep('Spleen', 2+2+2), rep('Spleen', 12), rep('Spleen', 12-4), rep(c('Spleen'), 8)),# rep('Thymus', 2), rep('Thymus', 7), 'Thymus', 
                                   c(rep('WT', 7+8-2), rep('WT', 12), 'Cre-', 'Cre-', 'Cre-', 'Cre-', 'WT', 'WT', 'WT', 'WT', rep('WT', 16-8)), # 'ilin28b', 1+ 'Cre+', 'Cre+', 'Cre+', 'Cre+', +7
                                   c(#rep('NEB', 2),
                                     rep('Illumina', 7),
                                     rep(c('Illumina', 'NEB'), 2), rep('Illumina', 4-2),
                                     rep('NEB', 12+12-4+16-8)),# +7
                                   c(#0.122, 0.071,
                                     0.475, 0.545, 0.526, 0.511, 0.446, 0.446, 0.469,
                                     0.647, 0.417, 0.481, 0.226, 0.629, 0.461, #0.692, 0.637, 
                                     0.6385, 0.6385, 0.6385, 0.4567, 0.4567, 0.4567, 0.708, 0.708, 0.708, 0.633, 0.633, 0.633,
                                     #0.608, 0.707, 0.67, 0.714, 0.594, 0.642, 0.66,
                                     0.376, 0.457, 0.559, 0.605, 0.477, 0.542, 0.548, 0.457, # , 0.323, 0.612, 0.495, 0.599
                                     0.264, 0.225, 0.49, 0.329, 0.333, 0.284, 0.39, 0.373), # mean of replicates (not including Clean_TN_rep1), # 0.593, 
                                   c(rep('Clean', 7+8-2+6), rep('Dirty', 6), rep('Clean', 12-4+16-8)) # 2+ rep('Clean', 4), rep('Dirty', 3), 
)
colnames(pheno) = c("Sample", "Source", "Age", "CellType", "Organ", "Treatment", "Kit", "MappingRate", "Clean") # 
rownames(pheno) = pheno$Sample
pheno$Age = factor(pheno$Age)
pheno = as.matrix(pheno)

# data from anamiR package for debug
# pheno2 = pheno[1:30,]
# rownames(pheno2) <- colnames(mrna)
# data(mrna)
# data(mirna)
# data(pheno.mirna)
# data(pheno.mrna)
# 
# mrna_se <- SummarizedExperiment(
#   assays = SimpleList(counts=mrna),
#   colData = pheno.mrna)
# 
# mirna_se <- SummarizedExperiment(
#   assays = SimpleList(counts=mirna),
#   colData = pheno.mirna)
# 

mrna_se <- SummarizedExperiment(
  assays = SimpleList(counts=as.matrix(mrna)),
  colData = as.matrix(pheno))

mirna_se <- SummarizedExperiment(
  assays = SimpleList(counts=as.matrix(mirna)),
  colData = as.matrix(pheno))



# function adapted for my input (original: https://rdrr.io/bioc/anamiR/src/R/differExp_discrete.R)
# previous:
# gp1 <- which(pheno_data == levels(pheno_data)[1])
# gp2 <- which(pheno_data == levels(pheno_data)[2])
# current:
# gp1 <- which(pheno_data == levels(factor(pheno_data))[1])
# gp2 <- which(pheno_data == levels(factor(pheno_data))[2])
# reason for this change: 
# pheno_data is generated during the function, which is the column in colData that contains the group information.
# levels(pheno_data) is supposed to get the names of the two groups, and then gp1 and gp2 is used to find samples that belong to each group. Works for the input of their package.
# but for my input, levels(pheno_data) is not able to get the names of the two groups (instead of Adult and Neo, just NULL), and gp1 and gp2 is not able to be used to find samples that belong to each group.
# I added factor function so levels(factor(pheno_data)) can get the names of the two groups (Adult and Neo). 
# I did 'factor' for this column of colData before calling the function (line 44 of this script), but it doesn't help.

differExp_discrete_edit <- function(
    se,
    class,
    method = c("t.test",
               "limma",
               "wilcox.test",
               "DESeq"),
    limma.trend = FALSE,
    t_test.var = FALSE,
    log2 = FALSE,
    p_value.cutoff = 0.05,
    p_adjust.method = "BH",
    logratio = 0.5
) {
  
  data <- SummarizedExperiment::assays(se)[[1]]
  
  if (log2 %in% "TRUE") {
    data <- log2(data)
  }
  
  method <- match.arg(method)
  pheno_data <- SummarizedExperiment::colData(se)[[class]]
  
  if (!is.null(pheno_data)) {
    pheno_data <- t(pheno_data)
    # seperate group
    if (method == "limma") {
      group <- t(pheno_data)
      type <- as.character(unique(unlist(group)))
      levels(group)[levels(group) == type[1]] <- 0
      levels(group)[levels(group) == type[2]] <- 1
    }
    gp1 <- which(pheno_data == levels(factor(pheno_data))[1]) # before edit: gp1 <- which(pheno_data == levels(pheno_data)[1])
    gp2 <- which(pheno_data == levels(factor(pheno_data))[2]) # before edit: gp2 <- which(pheno_data == levels(pheno_data)[2])
    
    #Fold Changes and mean
    if (length(gp1) == 1) {
      mean_gp1 <- data[, gp1]
    } else {
      mean_gp1 <- apply(data[, gp1], 1, mean)
    }
    if (length(gp2) == 1) {
      mean_gp2 <- data[, gp2]
    } else{
      mean_gp2 <- apply(data[, gp2], 1, mean)
    }
    
    if (method != "DESeq") {
      FC <- mean_gp1 - mean_gp2
      p_value <- vector(mode = "numeric", length = nrow(data))
    }
    
    # t.test
    if (method == "t.test") {
      t_test <- function (da, gp_1, gp_2) {
        stats::t.test(da[gp_1], da[gp_2], var.equal = t_test.var)[["p.value"]]
      }
      p_value <- apply(data, 1, t_test, gp1, gp2)
    }
    
    # wilcoxon
    if (method == "wilcox.test") {
      wilcoxon <- function (da, gp_1, gp_2) {
        stats::wilcox.test(da[gp_1], da[gp_2])[["p.value"]]
      }
      p_value <- apply(data, 1, wilcoxon, gp1, gp2)
    }
    
    # limma  (trend)
    if (method == "limma") {
      design <- stats::model.matrix(~ 0 + group)
      fit <- limma::lmFit(data, design)
      mc <- limma::makeContrasts("group0 - group1", levels = design)
      fit2 <- limma::contrasts.fit(fit, mc)
      eb <- limma::eBayes(fit2, trend = limma.trend)
      p_value <- eb[["p.value"]]
    }
    
    #DESeq
    if (method == "DESeq") {
      tmp <- as.formula(paste("~", class))
      dds <- DESeq2::DESeqDataSet(se, design = tmp)
      dds <- DESeq2::DESeq(dds)
      res <- DESeq2::results(dds)
      p_value <- res[["pvalue"]]
      p_adjust <- res[["padj"]]
      FC <- res[["log2FoldChange"]]
    }
    
    # output
    if (method != "DESeq") {
      p_adjust <- stats::p.adjust(p = p_value, method = p_adjust.method)
    }
    idx <- which(p_adjust < p_value.cutoff)
    DE_data <- data[idx, ]
    DE_data <- cbind(DE_data, FC[idx], p_value[idx], p_adjust[idx],
                     mean_gp1[idx], mean_gp2[idx])
    len_col <- ncol(DE_data)
    colnames(DE_data)[(len_col - 4):len_col] <- c("log-ratio",
                                                  "P-Value",
                                                  "P-adjust",
                                                  "mean_case",
                                                  "mean_control")
    FC_rows <- abs(DE_data[, len_col - 4])
    DE_data <- DE_data[which(FC_rows > logratio), ]
    return(DE_data)
  }
}

mrna_d <- differExp_discrete_edit(se = mrna_se,
                             class = "Age", method = "t.test",
                             t_test.var = FALSE, log2 = FALSE,
                             p_value.cutoff = 0.05,  logratio = 0.5
)

mirna_d <- differExp_discrete_edit(se = mirna_se,
                              class = "Age", method = "t.test",
                              t_test.var = FALSE, log2 = FALSE,
                              p_value.cutoff = 0.05,  logratio = 0.5
)

cor <- negative_cor(mrna_data = mrna_d, mirna_data = mirna_d,
                    method = "pearson", cut.off = -0.5)

# merge with miRNA names

mir_family_info <- read.table("targetscan/mouse/miR_Family_Info.txt", head = T, fill = T, sep = "\t") 
hsa_mir_family_info <- mir_family_info[substr(mir_family_info[,'MiRBase.ID'] ,1 ,3) == "mmu", ] 
hsa_con_mir_family_info <- hsa_mir_family_info[hsa_mir_family_info[,'Family.Conservation.'] >= 2, ] 
seed_con <- hsa_con_mir_family_info[,c('Seed.m8', 'MiRBase.ID')] 
mirna_seeds_con <- aggregate(seed_con[,2], by = list(seed = seed_con$Seed), FUN = paste) 
colnames(mirna_seeds_con) <- c("seed", "miRNA") 
mirna_seeds_con$miRNA_name = c('miR-106ab,17,20ab,93-5p,6383','miR-141,200a-3p','miR-132,212-3p','miR-451a','miR-490-3p','miR-191-5p','miR-124-3p.1','miR-18ab-5p','miR-291a,294,295,302abd-3p','miR-200bc,429-3p','miR-216a-5p','miR-216b-5p','miR-365-3p','miR-101a-3p.1','miR-144-3p','miR-181abcd-5p','miR-140-3p.2,miR-497b','miR-100,99ab-5p','miR-10ab-5p','miR-217-5p','miR-193ab-3p','miR-383-5p.2','miR-29abc-3p','miR-15ab,16,1907,195a,322,497a-5p,miR-195b,6342,6353,6419','miR-129-1,2-3p','miR-22-3p','miR-21a,590-5p,miR-21c','miR-196ab-5p','miR-130ab,301ab-3p,miR-130c,6341,6389,721','miR-302c-3p','miR-140-5p,miR-876-3p','miR-142a-5p','miR-425-5p,miR-489-3p','miR-183-5p','miR-135ab-5p','miR-455-5p,miR-5129-3p','miR-25,363,367,92ab-3p,miR-32-5p','miR-128-3p,miR-6539','miR-802-5p','miR-199ab-3p','miR-455-3p.1','miR-148ab,152-3p','miR-140-3p.1','miR-338-3p','miR-199ab-5p','miR-125ab,351-5p,miR-6367,6394','miR-205-5p','miR-212-5p','miR-551b-3p','miR-126a-3p.1','miR-187-3p','miR-139-5p','miR-150-5p,miR-5127','miR-9-5p','miR-203-3p.1','miR-146ab-5p','miR-143-3p','let-7abcdefgik,98-5p,miR-1961','miR-190ab-5p','miR-383-5p.1','miR-219a-5p','miR-103,107-3p','miR-214-5p','miR-221,222-3p,miR-1928','miR-138-5p','miR-7ab-5p','miR-1a,206-3p,miR-1957b,6349,6382','miR-184-3p','miR-122-5p','miR-31-5p','miR-34abc,449ac-5p,-miR-449b','miR-24-3p,miR-5124b,6361,6369,6410,6413','miR-193a-5p','miR-30abcde,384-5p','miR-194-5p','miR-126a-3p.2','miR-142a-3p.1','miR-223-3p','miR-19ab-3p','miR-208ab-3p','miR-499-5p','miR-124,5624-3p,miR-6540-5p','miR-155-5p','miR-101ab-3p','miR-142a-3p.2','miR-137-3p','miR-26ab-5p','miR-27ab-3p','miR-23ab-3p','miR-145a-5p,miR-145b','miR-204,211-5p,miR-7670-3p','miR-202-5p','miR-203-3p.2','miR-192,215-5p','miR-455-3p.2,miR-682','miR-153-3p','miR-33-5p','miR-183-5p.2','miR-133a-3p.1','miR-147-3p','miR-210-3p','miR-218-5p,miR-7002-3p','miR-182-5p','miR-96-5p','miR-133ab-3p,miR-133c','miR-375-3p','miR-129-5p')

mirna_seeds_con$miRNA = NULL

cor_miRNAname = merge(mirna_seeds_con, cor, by.x = 'seed', by.y = 'miRNA')
write.table(cor_miRNAname, 
            'make_figures/benchmark/anamiR/res.txt',
            sep = '\t',
            quote = F, row.names = F)
# mir29 targets
mir29targets = cor_miRNAname[cor_miRNAname$miRNA_name == 'miR-29abc-3p','Gene']
# no let7 targets found
