# DNAzyme Enrichment

Analysis scripts for DNAzyme enrichment NGS data processing and motif/position enrichment summaries.

This repository is organized for sharing code and workflow only.
Large raw data and generated outputs are intentionally excluded from Git tracking.

## Project Layout

- `script/`: main analysis scripts and notebooks
  - `bash_per.sh`: run one sample/round pipeline
  - `bash_grand.sh`: run all configured rounds
  - `MotifFreq_byRounds_v1.py`: motif-position counting and heatmap generation
  - `MotifFreq_Enrichment_v2.ipynb`: enrichment exploration over rounds notebook
- `data/`: input data
- `results/`: analysis outputs
- `log/`: run logs (ignored)
- `environment.yml`: reproducible conda environment

## Data Assumptions

The shell workflow expects this input structure:

- `data/aKG/aKG18_L001/...fastq.gz`
- `data/aKG/aKG19_L001/...fastq.gz`
- `data/aKG/aKG21_L001/...fastq.gz`
- `data/nKG/nKG19_L001/...fastq.gz`
- `data/nKG/nKG20_L001/...fastq.gz`
- `data/nKG/nKG22_L001/...fastq.gz`

Naming convention expected by `bash_per.sh`:
- paired files include `R1_001.fastq.gz` and `R2_001.fastq.gz`

Configuration Note: 
- Update the sequence variables in script/bash_per.sh if your DNAzyme construct or flanking adapters differ from the default design (Line # 81-89).

# Define DNAzyme sequence.

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

