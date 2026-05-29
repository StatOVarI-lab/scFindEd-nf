# scFindEd - nf

A Nextflow pipeline for preprocessing targeted amplicon sequencing data. 

The pipeline processes barcoded paired-end FASTQ files together with an amplicon design file to produce aligned reads and amplicon-level read summaries.

## Overview
1. Build reference
2. Create amplicon list
3. Extract read names
4. Align reads
5. Generate QC stats
6. Extract amplicon-level read information

## Requirements
* Nextflow (>22 recommended)
* Java (required for Nextflow)
* Conda

## Getting Started
### 1. Install Nextflow and Conda

See: https://www.nextflow.io/docs/latest/install.html 

### 2. Clone the repository
```
git clone https://github.com/StatOVarI-lab/scFindEd-nf.git
cd scFindEd-nf
```

## Running the pipeline

Example:

```
nextflow run main.nf -profile conda \
  --design tests/input/9824-design-summary.csv \
  --reads 'tests/input/*_R{1,2}_*.fastq' \
  --snp_info tests/input/20260310_Pilot1_SNP_final.xlsx \
  --outdir results \
  -resume
```
Note: Test datasets are not included in the repository due to file size.

## Key parameters
| Parameter   | Description                                   |
|-------------|-----------------------------------------------|
| `--design`  | Amplicon design CSV file                      |
| `--reads`   | FASTQ file pattern for paired reads; assumes FASTQ are barcoded          |
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

```
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
│   └── sample_ampliconID.tsv
│   └── sample_ampliconID.tsv
│   └── sample_ampliconID.tsv
│
├── readnames/
│   └── sample.readnames.tsv
│ 
├── design/
│   └── 20260310_Pilot1_SNP_final.with_AmpIDs.csv
│   └── 9824-design-summary.joined.csv
│   └── 6915_design_summary.joined.csv
│   └── 9824-design-summary.with_amplicon_seq.csv
│   └── 6915_design_summary.with_amplicon_seq.csv
│   └── old_pilot_SNPs.with_AmpIDs.csv
│
├── fragment_mutations/
│   └── sample_name/
│           └── amplID.fragment_mutations.tsv
│           └── amplID.fragment_mutations.tsv
│           └── amplID.fragment_mutations.tsv
│           └── amplID.fragment_mutations.tsv
│
├── qc/
│    └── sample.flagstat.txt
│    └── sample.idxstats.txt
│   
│
└── reference_amplicon/
    ├── reference.fa
    ├── reference.fa.fai
    └── reference.dict
```

### Mutation evidence output
The Nextflow workflow does not perform final mutation summarisation. Instead, it produced intermediate mutation-evidence files that can be used by downstream analysis workflows such as scFindEd R package. 

For each ampliconm, the workflow outputs key information required to summarise mutations (`results/fragment_mutations/sample/file.tsv`): 

| Column | Description | Downstream use | Example / cutoff |
|---|---|---|---|
| `sample_id` | Sample identifier | Separate biological samples | `MB_barcoded` |
| `cell_barcode` | Cell barcode sequence | Group mutations per cell | `AACGTT...-1` |
| `amp_id` | Amplicon identifier | Link mutation to target region | `AMPL904013` |
| `r1_pos` | R1 alignment start | Determine overlap coordinates | `5` = read starts at amplicon position 5 |
| `r1_mapq` | R1 mapping quality | Filter low-confidence alignments | Often keep `MAPQ >= 30` |
| `r1_cigar` | R1 alignment structure | Map read to reference | `96M`, `38S60M` |
| `r1_flag` | R1 SAM flag | Filter alignment types | `99` = proper paired R1 |
| `r1_nm` | R1 edit distance | Count total read differences | `NM = 0` reference match |
| `r1_md` | R1 mismatch string | Identify mismatch positions | `73T22` = mismatch after 73 bases |
| `r1_seq` | R1 read sequence | Extract observed bases | Used to obtain mutation base |
| `r2_pos` | R2 alignment start | Determine overlap coordinates | `81` |
| `r2_mapq` | R2 mapping quality | Filter low-confidence alignments | Often keep `MAPQ >= 30` |
| `r2_cigar` | R2 alignment structure | Map read to reference | `150M` |
| `r2_flag` | R2 SAM flag | Filter alignment types | `147` = proper paired R2 |
| `r2_nm` | R2 edit distance | Count total read differences | `NM > 0` indicates differences |
| `r2_md` | R2 mismatch string | Identify mismatch positions | `51^A99` = deletion |
| `r2_seq` | R2 read sequence | Extract observed bases | Used to obtain mutation base |
| `n_fragments` | Fragment support count | Weight mutation support | Often require `>= 2–5` fragments |

### Potential steps for processing mutations from `results/fragment_mutations/` directory:

The mutation-evidence step starts from grouped fragment patterns, where each row represents one unique R1/R2 alignment pattern and `n_fragments` records how many fragments share that pattern.

1. **Read grouped fragment patterns**  
   Import the grouped fragment-pattern table generated from the aligned BAM files.

2. **Filter alignments using SAM flags**  
   Decode the SAM `FLAG` field and retain only:
   - paired reads (`FLAG & 1`)
   - properly paired reads (`FLAG & 2`)
   - primary alignments
   
   Remove:
   - secondary alignments (`FLAG & 256`)
   - supplementary alignments (`FLAG & 2048`)
   
   Incomplete pairs (`R1_only` and `R2_only`) are separated by requiring both R1 and R2 to be present after pairing.

3. **Filter rows by NM values**  
   Separate rows where `r1_nm > 0` or `r2_nm > 0`, because rows with `NM = 0` for both reads have no mismatches to extract.

4. **Parse CIGAR strings**  
   The `CIGAR` string is split into operation lengths and operation types, such as `96M`, `38S60M`, or `51M1D99M`. This is used to track how read positions relate to amplicon-reference positions, including soft clipping, insertions and deletions.

5. **Parse MD strings**  
   The `MD` string is used to identify mismatch positions and reference bases. For example, `73T22` means 73 matching bases, then a mismatch where the reference base is `T`, followed by 22 matching bases.

6. **Map mismatch to amplicon position**  
   The mismatch position from `MD` is converted into an amplicon coordinate using the read alignment start position, for example `amplicon_pos = r1_pos + ref_index - 1`.

7. **Extract observed read base**  
   The `CIGAR` string is used to map the mismatch reference position back to the correct position in the read sequence. The observed base is then extracted from `r1_seq` or `r2_seq`.

8. **Create mutation call**  
   The mutation is reported as `ref_base > read_base`, for example `T>A`.

9. **Compare with mate read**  
   For each mutation detected in R1, the script checks what base R2 has at the same amplicon position, and vice versa.

10. **Classify paired-read support**  
   Each mutation is classified as:
   - `R1_nonoverlap_only`: mutation detected in R1, but R2 does not cover the position
   - `R2_nonoverlap_only`: mutation detected in R2, but R1 does not cover the position
   - `R1_R2_agree`: both reads cover the position and show the same mutated base
   - `R1_R2_disagree`: both reads cover the position but show different bases

11. **Propagate fragment support**  
   The `n_fragments` value is carried forward so downstream workflows can weight each mutation by the number of identical fragments supporting that evidence pattern.

## Pipeline Reports
The pipeline does NOT automatically generates execution reports. 

## Development Notes
The pipeline is implemented in Nextflow DSL2.

Intermediate files are stored in the `work/ `directory and final results are published to the directory specified by `--outdir`.
