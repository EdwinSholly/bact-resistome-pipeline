#!/bin/bash
# ==============================================================================
# BACTERIAL GENOMIC RESISTOME PROFILING PIPELINE
# Core Modular Execution Script
# ==============================================================================

# Exit immediately if a command exits with a non-zero status
set -e

echo "=== Starting Genonmic Resistome Pipeline ==="

# 1. Quality Control
mkdir -p qc_reports
fastqc raw_data/sample_R1.fastq.gz raw_data/sample_R2.fastq.gz -o qc_reports/

# 2. Adapter & Quality Trimming
mkdir -p trimmed
trimmomatic PE -threads 4 \
  raw_data/sample_R1.fastq.gz raw_data/sample_R2.fastq.gz \
  trimmed/sample_R1_paired.fq.gz trimmed/sample_R1_unpaired.fq.gz \
  trimmed/sample_R2_paired.fq.gz trimmed/sample_R2_unpaired.fq.gz \
  ILLUMINACLIP:$CONDA_PREFIX/share/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 \
  SLIDINGWINDOW:4:20 MINLEN:50

# 3. Reference Indexing & Alignment
bwa index reference.fna
mkdir -p alignments
bwa mem -t 4 -R "@RG\tID:run1\tSM:sample\tPL:ILLUMINA\tLB:lib1\tPU:lane1" \
  reference.fna trimmed/sample_R1_paired.fq.gz trimmed/sample_R2_paired.fq.gz \
  > alignments/sample.sam

# 4. Coordinate Sorting & Binary Compression
samtools sort -@ 4 alignments/sample.sam -o alignments/sample_sorted.bam
samtools index alignments/sample_sorted.bam

# 5. PCR Duplicate Removal & Flagging
picard MarkDuplicates \
  -I alignments/sample_sorted.bam -O alignments/sample_dedup.bam \
  -M alignments/dup_metrics.txt -REMOVE_DUPLICATES false
samtools index alignments/sample_dedup.bam

# 6. Reference Pre-processing for GATK
samtools faidx reference.fna
gatk CreateSequenceDictionary -R reference.fna -O reference.dict

# 7. De Novo Variant Calling (GATK HaplotypeCaller)
mkdir -p variants
gatk HaplotypeCaller -R reference.fna -I alignments/sample_dedup.bam -O variants/raw.vcf

# 8. Hard Quality Filtration Matrix
gatk VariantFiltration -R reference.fna -V variants/raw.vcf -O variants/tagged.vcf \
  --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0" --filter-name "BASIC_FILTER"
gatk SelectVariants -R reference.fna -V variants/tagged.vcf --exclude-filtered -O variants/clean.vcf

# 9. Functional Annotation via SnpEff
snpEff -Xmx4g -v Staphylococcus_aureus variants/clean.vcf > variants/annotated.vcf

# 10. Reference-Free Consensus Generation
bgzip -c variants/clean.vcf > variants/clean.vcf.gz
tabix -p vcf variants/clean.vcf.gz
bcftools consensus -f reference.fna variants/clean.vcf.gz > variants/whole_genome_consensus.fasta

# 11. Mobilome Recovery Net (Extracting Unmapped Reads via SAM Flag 4)
samtools fasta -f 4 alignments/sample_sorted.bam > variants/unmapped_salvage.fasta

echo "=== Pipeline Workflow Execution Completed Successfully ==="
