### samll RNAseq pipeline

# quality check
fastqc -t 12 -f fastq -o qc fastq/*

# adapter trimming
mapper.pl config.txt -d -e -h -k AGATCGGAAGAGCACACGTCT -l 18 -m -s reads.fa -v

# mapping to miRNA sequences
quantifier.pl -p /path/to/mirbase/22/hairpin_ref_mmu.fa -m /path/to/mirbase/22/mature_ref_mmu.fa -r reads.fa -t mmu -d -e 0 -W 2> quantifier-W.log



### RNAseq pipeline

# quality check
fastqc -t 12 -f fastq -o qc fastq/*.fastq.gz

# cut adapters
trim_galore --fastqc_args "--outdir qc_cutadapt" --output_dir trim_galore --cores 7 fastq/*.fastq.gz 2> trim_galore/trim_galore.log

# mapping
for f in `ls -1 *_trimmed.fq.gz | sed 's/_trimmed.fq.gz//' `
do
hisat2 -x /path/to/mouse/mm10/genome -U ${f}_trimmed.fq.gz -S ../hisat2/${f}.sam -p 12 2> ../hisat2/${f}_summary.txt
done

# count reads
featureCounts -a /path/to/mouse/gencode_annotation/gencode.vM21.annotation.names.saf -F 'SAF' -s 0 -Q 50 -p -o counts.txt hisat2/*.sam 2> featureCounts.log

