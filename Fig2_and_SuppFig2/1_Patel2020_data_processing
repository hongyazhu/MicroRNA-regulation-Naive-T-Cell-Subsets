# Data used for this fig is from Patel et al., 2020 
# RNAseq: GSE142895

fastqc -t 12 -f fastq -o qc fastq_data/*

for i in $(ls *.fastq | sed 's/.fastq//'); do fastx_clipper -i ${i}.fastq -l 15 -Q 33 -a GATCGGAAGAGCACACGTCTGAACTCCAGTC -v -o ../trim/${i}\_trimmed.fastq > ../trim/${i}_trim.log & done

for f in `ls -1 *_trimmed.fastq.gz | sed 's/_trimmed.fastq.gz//' `
do
hisat2 -x /path/to/human/hg38/genome -U ${f}_trimmed.fastq.gz -S ../hisat2/${f}.sam -p 10 2> ../hisat2/${f}_summary.txt
done

# count exon and intron reads*
featureCounts -T 10 -a /path/to/human/gencode_annotation/gencode.v30.annotation.gene.names.noOverlap.saf -F 'SAF' -s 0 -Q 50 -o carpRNAseq_counts_gb.txt hisat2/*.sam 2> carpRNAseq_featureCounts_gb.log
featureCounts -T 10 -a /path/to/human/gencode_annotation/gencode.v30.annotation.names.noOverlapGenes.saf -F 'SAF' -s 0 -Q 50 -o carpRNAseq_counts_exon.txt hisat2/*.sam 2> carpRNAseq_featureCounts_exon.log
* saf files generated as below


# Generating noOverlap saf files
# from gtf to saf (with gene names)
# select lincRNA/miRNA/protein_coding type of genes/transcripts
awk '$3 == "exon"' gencode.v30.annotation.gtf | cut -f 1,4,5,6,7,9 | grep "lincRNA\|miRNA\|protein_coding" | awk --field-separator=";" '$3 ~ /protein_coding|lincRNA|m
iRNA/ ' | awk --field-separator=";" '$5 ~ /protein_coding|lincRNA|miRNA/ ' > gencode.v30.annotation.tmp.gtf
# get names, dots, strand to form .bed file
cut -f 6 gencode.v30.annotation.tmp.gtf | awk 'FS="\"; "{print $4}' | sed 's/gene_name \"//' > gene_names.txt
nano gene_names.txt # set first row to MIR6859-1
cut -f 4 gencode.v30.annotation.tmp.gtf > dots.txt
cut -f 5 gencode.v30.annotation.tmp.gtf > strand.txt
cut -f 1-3 gencode.v30.annotation.tmp.gtf | paste -d "\t" - gene_names.txt | paste -d "\t" - dots.txt | paste -d "\t" - strand.txt > gencode.v30.annotation.names.bed
# sort and merge the .bed file to get saf file
sort -k1,1 -k2,2n gencode.v30.annotation.names.bed > gencode.v30.annotation.names.sort.bed
bedtools merge -i gencode.v30.annotation.names.sort.bed -s -c 4,6 -o distinct | awk '!/,/' | awk 'BEGIN {FS="\t";OFS="\t"} {print $4,$1,$2,$3,$5}' > gencode.v30.annotation.names.saf

# for genes (including introns) from gtf to saf (with gene names)
awk '$3 == "gene"' gencode.v30.annotation.gtf | cut -f 1,4,5,6,7,9 | grep "lincRNA\|miRNA\|protein_coding" | awk --field-separator=";" '$2 ~ /protein_coding|lincRNA|miRNA/ ' > gencode.v30.annotation.gene.tmp.gtf
cut -f 6 gencode.v30.annotation.gene.tmp.gtf | awk 'FS="\"; "{print $3}' | sed 's/gene_name \"//' > gene.gene_names.txt
vim gene.gene_names.txt #set first row to MIR6859-1
cut -f 4 gencode.v30.annotation.gene.tmp.gtf > gene.dots.txt
cut -f 5 gencode.v30.annotation.gene.tmp.gtf > gene.strand.txt
cut -f 1-3 gencode.v30.annotation.gene.tmp.gtf | paste -d "\t" - gene.gene_names.txt | paste -d "\t" - gene.dots.txt | paste -d "\t" - gene.strand.txt > gencode.v30.annotation.gene.names.bed
sort -k1,1 -k2,2n gencode.v30.annotation.gene.names.bed > gencode.v30.annotation.gene.names.sort.bed
# not running: bedtools merge -i gencode.vM21.annotation.gene.names.sort.bed -s -c 4,6 -o distinct | awk '!/,/' | awk 'BEGIN {FS="\t";OFS="\t"} {print $4,$1,$2,$3,$5}' > gencode.vM21.annotation.gene.names.saf

# get saf with genes that do not overlap with each other (help from https://stackoverflow.com/questions/43432149/filter-overlapping-entries-in-bed-file):
# gene body:
# 1) Count the overlaps in the original input
bedtools merge -i gencode.v30.annotation.gene.names.sort.bed -c 1 -o count > gencode.v30.annotation.gene.names.sort.bed.counted
# 2) Filter out only those rows that do not overlap with anything
awk '/\t1$/{print}' gencode.v30.annotation.gene.names.sort.bed.counted > gencode.v30.annotation.gene.names.sort.bed.filtered
# 3) Intersect it with the original input and keep only those original rows that were found after filtering as well
bedtools intersect -a gencode.v30.annotation.gene.names.sort.bed -b gencode.v30.annotation.gene.names.sort.bed.filtered -wa > gencode.v30.annotation.gene.names.sort.noOverlap.bed
bedtools merge -i gencode.v30.annotation.gene.names.sort.noOverlap.bed -s -c 4,6 -o distinct | awk '!/,/' | awk 'BEGIN {FS="\t";OFS="\t"} {print $4,$1,$2,$3,$5}' > gencode.v30.annotation.gene.names.noOverlap.saf

# exon:
# Intersect filtered list from gene file with the exon bed file and keep only those rows that were found after filtering as well (to avoid removing genes whose exons overlap with each over)
bedtools intersect -a gencode.v30.annotation.names.sort.bed -b gencode.v30.annotation.gene.names.sort.bed.filtered -wa > gencode.v30.annotation.names.sort.noOverlapGenes.bed
bedtools merge -i gencode.v30.annotation.names.sort.noOverlapGenes.bed -s -c 4,6 -o distinct | awk '!/,/' | awk 'BEGIN {FS="\t";OFS="\t"} {print $4,$1,$2,$3,$5}' > gencode.v30.annotation.names.noOverlapGenes.saf
