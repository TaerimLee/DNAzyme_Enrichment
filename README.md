# DNAzyme Enrichment

Analysis scripts for DNAzyme enrichment NGS data processing and motif/position enrichment summaries.

This repository is organized for sharing code, workflow, and results only.
Large raw sequencing data and outputs are intentionally excluded from Git tracking.

## Project Layout

- `script/`: main analysis scripts and notebooks
  - `bash_per.sh`: run one sample/round pipeline
  - `bash_grand.sh`: run all configured rounds
  - `MotifFreq_byRounds_v1.py`: motif-position counting and heatmap generation
  - `MotifFreq_Enrichment_v2.ipynb`: enrichment exploration over rounds notebook
- `data/`: input data
- `results/`: analysis outputs
- `log/`: run logs
- `environment.yml`: reproducible conda environment

## Input Data

The shell workflow expects this input structure:
- aKG
  - `data/aKG/aKG18_L001/aKG18_S6_L001_R1_001.fastq.gz`
  - `data/aKG/aKG18_L001/aKG18_S6_L001_R2_001.fastq.gz`
  - `data/aKG/aKG19_L001/aKG19_S7_L001_R1_001.fastq.gz`
  - `data/aKG/aKG19_L001/aKG19_S7_L001_R2_001.fastq.gz`
  - `data/aKG/aKG21_L001/aKG21_S8_L001_R1_001.fastq.gz`
  - `data/aKG/aKG21_L001/aKG21_S8_L001_R2_001.fastq.gz`
- nKG
  - `data/nKG/nKG19_L001/nKG19_S9_L001_R1_001.fastq.gz`
  - `data/nKG/nKG19_L001/nKG19_S9_L001_R2_001.fastq.gz`
  - `data/nKG/nKG20_L001/nKG20_S10_L001_R1_001.fastq.gz`
  - `data/nKG/nKG20_L001/nKG20_S10_L001_R2_001.fastq.gz`
  - `data/nKG/nKG22_L001/nKG22_S11_L001_R1_001.fastq.gz`
  - `data/nKG/nKG22_L001/nKG22_S11_L001_R2_001.fastq.gz`

Naming convention expected by `bash_per.sh`:
- paired files include `R1_001.fastq.gz` and `R2_001.fastq.gz`


# Define DNAzyme sequence.
- Update the sequence variables in script/bash_per.sh if your DNAzyme construct or flanking adapters differ from the default design (Line # 81-89).
```
# 'Enzyme' sequence: 5'- TAATACGACTCACTATAGGG - N40 - ATCTGACGGTAACGCTATAGTGTCACCTAAATAGC -3’
#                        ***TACGACTCACTATAGGG - N40 - ATCTGACGGTAACGCTATAGTGTCACCTAAAT***
## Remove first (5`) and last (3`) 3 bases for buffering. Since cutadapt can remove a few more/less bases, we will allow 3 bases buffering.
## Define flanking sequence.

seq_5="TACGACTCACTATAGGG"
seq_3="ATCTGACGGTAACGCTATAGTGTCACCTAAAT"
seq_3_rc="ATTTAGGTGACACTATAGCGTTACCGTCAGAT" # Reverse complemented
```

## Quick Start

```bash
# From repository root
conda env create -f environment.yml
conda activate DNAzymeEnrichment

# Run one sample/round
bash script/bash_per.sh aKG 19

# Run all predefined rounds
bash script/bash_grand.sh

# After running all sample/round, run /script/MotifFreq_Enrichment_v2.ipynb for enrichment analysis.
```

