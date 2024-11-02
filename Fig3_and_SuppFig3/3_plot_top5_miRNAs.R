# plot top 5 miRNAs with most predicted targets: Fig 3B-C.
# make two plots: one contains miR-Inf targets; one contains miR-Inf targets that are also conserved.

net_ct02_4nets = read.table('analyze_network/ms2/bias5/prioritizeTargets/pos_regs_freqTFmRNAmorethan4.tsv',
                              header = T)
length(table(net_ct02_4nets$miRNA_name))
length(table(net_ct02_4nets$target))

outputpath = "make_figures/ct02_ms2_4ormore"

# no. of regulations for each miRNA (plotting mirnas with more than 50 predictions)
num_regs_permirna = as.data.frame(table(net_ct02_4nets$miRNA_name))
colnames(num_regs_permirna) = c('miRNA', 'No.regulations')
num_regs_permirna <- within(num_regs_permirna, 
                            miRNA <- factor(miRNA, levels=names(sort(table(net_ct02_4nets$miRNA_name), decreasing = F))))
num_regs_permirna_sort <- num_regs_permirna[order(num_regs_permirna$No.regulations, decreasing = TRUE),]
library(ggplot2)
ggplot(data=num_regs_permirna_sort[1:5,], aes(x=miRNA, y=No.regulations)) +
  geom_bar(stat="identity", fill="bisque3", width = 0.7)+
  theme_bw() + 
  coord_flip() + ylab('No. regulations') + xlab('') +
  ggtitle('Inferelator-predicted\n(smRNAseq4/5)')+
  theme(plot.title = element_text(hjust = 0.5))
ggsave(paste0(outputpath, "/num_regs_permirna_inf_smRNAseq4ormore.pdf"), width = 15, height = 8, units = "cm")

# no. of regulations for each miRNA, conserved regulations only
net_ct02_4nets_con = net_ct02_4nets[net_ct02_4nets$ts == 'TRUE',]
num_regs_permirna_con = as.data.frame(table(net_ct02_4nets_con$miRNA_name))
colnames(num_regs_permirna_con) = c('miRNA', 'No.regulations')
num_regs_permirna_con <- within(num_regs_permirna_con, 
                                miRNA <- factor(miRNA, levels=names(sort(table(net_ct02_4nets_con$miRNA_name), decreasing = F))))
num_regs_permirna_con_sort <- num_regs_permirna_con[order(num_regs_permirna_con$No.regulations, decreasing = TRUE),]
library(ggplot2)
ggplot(data=num_regs_permirna_con_sort[1:5,], aes(x=miRNA, y=No.regulations)) +
  geom_bar(stat="identity", fill="bisque3", width = 0.7)+
  theme_bw() + 
  coord_flip() + ylab('No. regulations') + xlab('') +
  ggtitle('Inferelator-predicted, conserved\n(smRNAseq4/5)')+
  theme(plot.title = element_text(hjust = 0.5))
ggsave(paste0(outputpath, "/num_regs_permirna_infcon_smRNAseqmorethan4.pdf"), width = 15, height = 8, units = "cm")
