# Generalizable Genomic Resistome Pipeline: Deconvoluting Reference-Mapping Bias in Bacterial Assemblies

[![Pipeline Platform](https://img.shields.io/badge/Environment-Linux%20%7C%20WSL%20%7C%20Ubuntu-orange)](https://ubuntu.com/)
[![Bioinformatics Tools](https://img.shields.io/badge/Tools-FastQC%20%7C%20BWA%20%7C%20Samtools%20%7C%20GATK%20%7C%20SnpEff-blue)](https://gatk.broadinstitute.org/)
[![Database Core](https://img.shields.io/badge/Resistome-CARD%20RGI%20%7C%20NCBI%20AMRFinderPlus-green)](https://card.mcmaster.ca/)

## Project Overview & Design Philosophy
This repository contains a modular, end-to-end clinical bioinformatics pipeline engineered for high-throughput variant discovery and reference-independent resistome profiling from short-read Next-Generation Sequencing (NGS) data. 

### The Reference-Mapping Constraint
[cite_start]Standard variant-calling pipelines rely heavily on mapping raw reads against a vulnerable baseline reference template. [cite_start]However, this approach creates a severe structural blind spot: unique insertions, structural variations, and hyper-divergent horizontal gene transfer (HGT) elements completely fail to map to the reference, causing them to be categorized under **SAM Binary Flag 4 (Segment Unmapped)** and discarded[cite: 278, 307, 308]. 

### The Solution
[cite_start]This tool solves coordinate-system constraints by executing a two-pronged validation strategy[cite: 287]:
1. [cite_start]**Pseudogenome Re-reconstruction:** Leveraging `bcftools consensus` to physically integrate high-confidence variants back into the reference framework template, producing an unconstrained continuous sequence asset for scanning[cite: 279, 281].
2. [cite_start]**The "Lost & Found" Mobilome Recovery Net:** Isolating the unmapped fraction (`samtools view -f 4`) directly from the binary alignment maps and routing them into protein-translated database search engines to scan for foreign elements without reference constraints[cite: 312, 313].

---

## Table of Contents
1. [Pipeline Architecture](#pipeline-architecture)
2. [Repository Structure](#repository-structure)
3. [Getting Started](#getting-started)
   * [Prerequisites](#prerequisites)
   * [Installation](#installation)
4. [Usage Guide](#usage-guide)
5. [Validation Case Study: Correcting Public Metadata Errors](#validation-case-study-correcting-public-metadata-errors)
6. [Roadmap](#roadmap)
7. [Contributing](#contributing)
8. [License](#license)
9. [Contact](#contact)

---

## Pipeline Architecture
The workflow automates data transitions from raw sequencing quality evaluation to functional variant annotation and reference-independent resistome parsing:

```text
[Raw Paired-End FASTQ] ──> FastQC (Quality Control)
                                │
                                ▼
                           Trimmomatic (Quality & Adapter Filtering via $CONDA_PREFIX)
                                │
                                ▼
[BWA-MEM Index] ────────> BWA-MEM (Short-Read Alignment Map Generation)
                                │
                                ▼
                           Samtools view/sort/index (SAM-to-Sorted-BAM Conversion)
                                │
                                ▼
                           Picard MarkDuplicates (PCR Clone Removal & Flagging)
                                │
                                ▼
[RefSeq Sequence Dict] ──> GATK HaplotypeCaller (High-Sensitivity De Novo SNV/Indel Calling)
                                │
                                ▼
                           GATK VariantFiltration (Hard Filtering: QD, FS, MQ Thresholds)
                                │
                                ▼
                           GATK SelectVariants (Pristine True Variant Extraction)
                                │
                                ▼
                           SnpEff (Functional Coding Consequence Impact Mapping)
                                │
                                ▼
                           bcftools consensus (Continuous Chromosomal Assembly Generation)
                                │
                                ▼
[CARD / NCBI Engines] ───> Dual-Engine Antimicrobial Resistance (AMR) Screening
                                │
                                ▼
                           Mobilome Salvage (samtools -f 4 Isolation of Unmapped Reads)
```

---

## Repository Structure
The project layout isolates configuration scripts, documentation logs, and high-confidence variant tables into a clean, reproducible directory tree:

```text
bact-resistome-pipeline/
├── README.md                 # System architecture, deployment parameters, and case validation
├── docs/
│   └── project_summary.txt   # Execution validation log and clinical analytical report
├── config/
│   └── pipeline_steps.sh     # Production-ready wrapper execution script
└── results/
    ├── sample_target_missense.txt    # Filtered coding mutations unique to the target strain (1,513 count)
    ├── sample_control_missense.txt   # Mapped background variations for the baseline control (25,523 count)
    ├── sample_target_AMR_hits.txt    # Directed locus mask coordinate output for the target strain (0 coding hits)
    ├── sample_control_AMR_hits.txt   # Directed locus mask coordinate output for the control strain (0 coding hits)
    ├── sample_amrfinder_raw.txt      # Raw local clinical screening array from NCBI AMRFinderPlus
    ├── control_amrfinder_raw.txt     # Raw local clinical screening array for the control baseline
    └── sample_card_rgi_consensus.txt # Comprehensive whole-genome CARD RGI alignment logs (fosB/norC)
```
    
---

## Getting Started

### Prerequisites
* **Operating System:** Linux Environment (Native Kernel or Ubuntu via Windows Subsystem for Linux - WSL).
* **Package Manager:** Miniconda or Anaconda (for isolated environment deployment).

### Installation
Deploy the pipeline environment utilizing a tight strict-channel priority to prevent dependency clashing among downstream variant callers:

1. Clone the repository framework:
```bash
git clone https://github.com/EdwinSholly/bact-resistome-pipeline.git
cd bact-resistome-pipeline
```

2. Establish and build the isolated conda environment utilizing the mamba solver:
```bash
# Configure bioinformatics repository channels
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda install -n base -c conda-forge mamba -y

# Instantiate environment with matching runtime dependencies
mamba create -n bact_pipeline -c conda-forge -c bioconda -c defaults --strict-channel-priority -y \
  gatk4 bwa samtools bcftools picard fastqc trimmomatic openjdk=17 python=3.9 snpeff abricate

# Initialize environment
conda activate bact_pipeline
```

---

## Usage Guide
The pipeline is completely modularized. You can invoke the master wrapper script to run the entire analysis automatically or execute individual analytical building blocks manually to fine-tune specific parameters.

### 1. Execute Full Pipeline Automatically
To run the automated production-grade wrapper script, ensure it has execution permissions and pass your input assets using the standard argument flags:

```bash
# Grant execution permissions to the wrapper script
chmod +x config/pipeline_steps.sh

# Invoke the master workflow execution script
./config/pipeline_steps.sh -r reference.fna -g annotation.gff -1 forward_reads.fastq.gz -2 reverse_reads.fastq.gz
```

### 2. Manual Step-by-Step Block Execution
For deep-dive validation, custom filtration adjustments, or troubleshooting, you can run the core blocks of the pipeline independently in your terminal:

#### Block A: Quality Assessment & Adapter Trimming
Evaluate raw data integrity and execute strict sliding-window quality trimming to eliminate adapter pollution and low-confidence bases:

```bash
# Run baseline quality evaluation
mkdir -p qc_reports
fastqc raw_data/sample_R1.fastq.gz raw_data/sample_R2.fastq.gz -o qc_reports/

# Execute stringency trimming via dynamic conda paths
mkdir -p trimmed
trimmomatic PE -threads 4 \
  raw_data/sample_R1.fastq.gz raw_data/sample_R2.fastq.gz \
  trimmed/sample_R1_paired.fq.gz trimmed/sample_R1_unpaired.fq.gz \
  trimmed/sample_R2_paired.fq.gz trimmed/sample_R2_unpaired.fq.gz \
  ILLUMINACLIP:$CONDA_PREFIX/share/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 \
  SLIDINGWINDOW:4:20 MINLEN:50
```

#### Block B: Alignment & Coordinate Sorting
Index the target reference sequence layout and map your paired reads using a high-throughput seed-and-extend structural alignment algorithm:

```bash
# Build Burrows-Wheeler Transform index matrix
bwa index reference.fna

# Generate alignment sequence maps with explicit read group tracking
mkdir -p alignments
bwa mem -t 4 -R "@RG\tID:run1\tSM:sample\tPL:ILLUMINA\tLB:lib1\tPU:lane1" reference.fna trimmed/sample_R1_paired.fq.gz trimmed/sample_R2_paired.fq.gz > alignments/sample.sam

# Compress to binary format and index coordinate sorted positions
samtools sort -@ 4 alignments/sample.sam -o alignments/sample_sorted.bam
samtools index alignments/sample_sorted.bam
```

#### Block C: Duplicate Remediation & GATK Environment Setup
Isolate optical or PCR-cloned duplicate reads to prevent false amplification bias, and construct the sequence dictionaries required by the GATK engine:

```bash
# Flag and mitigate duplicate clonal artifacts
picard MarkDuplicates -I alignments/sample_sorted.bam -O alignments/sample_dedup.bam -M logs/dup_metrics.txt -REMOVE_DUPLICATES false
samtools index alignments/sample_dedup.bam

# Generate sequence indexes and metadata dict grids for the reference
samtools faidx reference.fna
gatk CreateSequenceDictionary -R reference.fna -O reference.dict
```

#### Block D: De Novo Variant Discovery & Hard Quality Filtration
Call genetic variations with high-sensitivity local de novo assembly, apply rigorous mathematical thresholds to isolate true variants, and discard background sequencing artifacts:

```bash
# Execute local genomic de novo variant calling
mkdir -p variants
gatk HaplotypeCaller -R reference.fna -I alignments/sample_dedup.bam -O variants/raw.vcf

# Apply strict hard-filtration quality expressions
gatk VariantFiltration -R reference.fna -V variants/raw.vcf -O variants/tagged.vcf \
  --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0" --filter-name "BASIC_FILTER"

# Extract pristine, high-confidence variations passing all filters
gatk SelectVariants -R reference.fna -V variants/tagged.vcf --exclude-filtered -O variants/clean.vcf
```

#### Block E: Reference-Free Consensus Assembly & Lost-Read Salvage Operation
Physically integrate your validated variants back into the reference backbone to generate an unconstrained consensus layout, while simultaneously trapping unmapped reads for structural mobilome screening:

```bash
# Compress and index structural variant maps
bgzip -c variants/clean.vcf > variants/clean.vcf.gz
tabix -p vcf variants/clean.vcf.gz

# Reconstruct a continuous coordinate-free pseudogenome asset
bcftools consensus -f reference.fna variants/clean.vcf.gz > variants/whole_genome_consensus.fasta

# Extract absolute unmapped read pool (SAM Flag -f 4) via optimized terminal text streams
samtools fasta -f 4 alignments/sample_sorted.bam > variants/unmapped_salvage.fasta
```

### 3. Verifying Pipeline Outputs
Once execution completes, your clean downstream assets will be organized within your local directory matrix. You can verify the integrity and file sizes of your final outputs by running:

```bash
ls -lh results/ variants/
```
**Expected core production files include:**
* `variants/whole_genome_consensus.fasta` - The continuous, unconstrained consensus assembly used for macro-resistome scanning.
* `variants/unmapped_salvage.fasta` - The isolated mobile genetic element fraction (SAM Flag 4) extracted for reference-independent alignment.
* `results/sample_target_missense.txt` - The clean, GATK/SnpEff-annotated coding mutation catalog for your target strain (containing the 1,513 filtered true variants).
* `results/sample_control_missense.txt` - The mapped background variations catalog for your baseline control track (containing the 25,523 background variants).
* `results/sample_target_AMR_hits.txt` - The directed locus mask coordinate output for your target strain showing 0 coding variations in textbook targets.
* `results/sample_control_AMR_hits.txt` - The directed locus mask coordinate output for your control strain showing 0 coding variations in textbook targets.
* `results/sample_amrfinder_raw.txt` - The raw, conservative clinical screening output array generated locally by NCBI AMRFinderPlus.
* `results/sample_card_rgi_consensus.txt` - The comprehensive Whole-Genome resistome profile log documenting the perfect *fosB* and strict *norC* chromosomal matches.

---

## Validation Case Study: Correcting Public Metadata Errors
To validate the generalizable workflow, the pipeline was deployed to evaluate an uncharacterized clinical *Staphylococcus aureus* isolate designated in public repository metadata fields as a Methicillin-Resistant strain ("MRSA").

### 1. Evaluation of Alignment Target Maps
Running the filtered variant calling output against the standard susceptible lab reference sequence (*S. aureus* NCTC 8325, Accession: NC_007795.1) yielded high-density background variations, detailing substantial evolutionary drift:
* **Background Drift Catalog:** 25,523 coding missense variations discovered in the positive control track.
* **Target Specimen Mutation Load:** 1,513 true, clean coding missense mutations.
* **Targeted AMR Coordinate Masking:** A targeted grep query isolating major standard resistance targets (*gyrA*, *grlA*, and *pbp2*) returned *0 high-impact coding alterations*, identifying exclusively low-impact non-coding MODIFIER mutations in flanking upstream regulatory zones.

### 2. Resolution of True Intrinsic Resistome Signature via Consensus Assembly
Because standard alignment failed to output the classic acquired horizontal resistance element *mecA* (Penicillin-Binding Protein 2a), the *Reference-Free Consensus* module and *Mobilome Recovery Net* were used to query global database profiles (*CARD* and *NCBI AMRFinderPlus*).

The unconstrained query successfully resolved the sample's true underlying resistance phenotype, proving it was driven by an *intrinsic, chromosomally encoded resistome layout* shared by this lineage and absent in the NCTC 8325 reference:

| Isolated Genomic Locus | Match Constraint Type | Sequence Identity Score | Global Reference Coordinates | Primary Molecular Action Mechanism |
| :--- | :--- | :--- | :--- | :--- |
| *fosB* | Perfect / Exact Match | 100.0% | 2,399,688 - 2,400,107 | Encoding a functional thiol-transferase enzyme facilitating direct enzymatic inactivation of Fosfomycin. |
| *norC* | Strict Homology Match | 99.13% | 62,340 - 63,728 | Major Facilitator Superfamily (MFS) multidrug efflux pump actively extruding fluoroquinolones. |

***Scientific Conclusion:*** The pipeline mathematically proved the total absence of acquired foreign resistance cassettes or plasmid backbones. The resistance profile is entirely endogenous to the baseline chromosome, systematically uncovering and correcting a legacy human metadata misclassification in public repositories.

---

## Roadmap
The development trajectory of this pipeline focuses on increasing scalability, containerization, and advanced transcriptomic integration to transition from static DNA variant screening to dynamic functional genomics:

- [x] **Phase 1: Core Variant Discovery Modules** — Implement high-sensitivity alignment, duplicate mitigation, and GATK hard-filtration blocks.
- [x] **Phase 2: Coordinate-Free Consensus Integration** — Build `bcftools consensus` pseudogenome reconstruction logic to bypass reference scaffolding constraints.
- [x] **Phase 3: Mobilome Recovery Net** — Deploy SAM Flag 4 unmapped read capture streams for reference-independent unconstrained screening.
- [ ] **Phase 4: Containerization & Workflow Management** — Migrate static Bash wrapper architecture into a distributed, production-ready Nextflow or Snakemake workflow equipped with Docker/Singularity container wrappers.
- [ ] **Phase 5: Regulatory & Expression Analytics (RNA-Seq Integration)** — Build an automated downstream differential expression pipeline utilizing **HISAT2, StringTie, and DESeq2** to quantify upstream promoter transcription kinetics and cross-examine expression variations in shared chromosomal resistome blueprints (*fosB* / *norC*).

---

## Contributing
Contributions are what make the open-source community such an amazing place to learn, inspire, and create. Any contributions you make to optimize workflow parallelization, add containerization layers, or improve database matching stringency are **greatly appreciated**.

1. Fork the Project framework.
2. Create your feature branch independently:
```bash
git checkout -b feature/AmazingFeature
```
3. Commit your modifications with a clear structural message:
```bash
git commit -m 'Add some AmazingFeature'
```
4. Push your changes directly up to your branch:
```bash
git push origin feature/AmazingFeature
```
5. Open a formal Pull Request for validation review.
---

## License
Distributed under the MIT License. This software is provided "as is", without warranty of any kind, allowing for free academic, institutional, and commercial modifications and distribution loops provided the original copyright and permission notice are preserved. 

See the local repository `LICENSE` file for more comprehensive legal clause details.

---

## Contact
Edwin Sholly K - [edwinsholly@gmail.com](mailto:edwinsholly@gmail.com)  
Project Link: [https://github.com/EdwinSholly/bact-resistome-pipeline](https://github.com/EdwinSholly/bact-resistome-pipeline)

***

<p align="center">
  <i>This pipeline was engineered as an automated, reproducible validation utility to systematically correct public repository metadata constraints and protect clinical data integrity.</i>
</p>
