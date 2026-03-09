# scFindEd - nf

A Nextflow pipeline for preprocessing targeted amplicon sequencing data. 

The pipeline processes paired-end FASTQ files together with an amplicon design file to produce aligned reads and amplicon-level read summaries.

## Overview
1. Build reference
2. Create amplicon list
3. Extract read names and barcodes
4. Align reads
5. Generate QC stats
6. Extract amplicon-level read information

## Requirements
Nextflow (>22 recommended)
Java (required for Nextflow)
Conda

## Getting Started
### 1. Install Nextflow and Conda

Use: https://www.nextflow.io/docs/latest/install.html 

### 2. Clone the repository
```
git clone https://github.com/StatOVarI-lab/scFindEd-nf.git
cd scFindEd-nf
```

## Running the pipeline

Example:

```
nextflow run main.nf -profile conda \
  --design tests/input/6915_design_summary_extended_information.csv \
  --reads "tests/input/*_R{1,2}_*.fastq" \
  --outdir results \
  --threads 10 \
  --compress false
```

## Key parameters
| Parameter   | Description                                   |
|-------------|-----------------------------------------------|
| `--design`  | Amplicon design CSV file                      |
| `--reads`   | FASTQ file pattern for paired reads           |
| `--outdir`  | Directory where results will be written       |
| `--threads` | Number of threads available for tools         |
| `--compress`| Whether TSV outputs should be gzipped         |

## Required inputs

### Amplicon design CSV

The design file must contain at least the following columns:

| Column              | Description                               |
|---------------------|-------------------------------------------|
| `AmpID`             | Unique identifier for each amplicon      |
| `Amplicon sequence` | Full amplicon sequence                   |

Rows missing these fields are filtered out during processing.

### Paired FASTQ files

Input reads should follow a paired-end naming pattern such as:
`sample_R1_001.fastq`
`sample_R2_001.fastq`

The pipeline expects a glob pattern such as:
`data/*_R{1,2}_*.fastq`

## Output
Results are written to the directory specified by `--outdir`.

Example structure:

results/
│
├── align/
│   └── sample.bam
│   └── sample.bam.bai
│
├── amplicons/
│   └── amplicon_list.txt
│
├── amplicon_tsv/
│   └── sample_ampliconID.tsv
│
├── readnames/
│   └── sample.readnames.tsv
│
├── qc/
│   
│
└── reference_amplicon/
    ├── reference.fa
    ├── reference.fa.fai
    └── reference.dict

## Pipeline Reports
The pipeline does NOT automatically generates execution reports. 

## Development Notes
The pipeline is implemented in Nextflow DSL2 and modularised into reusable process modules.

Intermediate files are stored in the `work/ `directory and final results are published to the directory specified by `--outdir`.
