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
  --design /path/to/design_summary.csv \
  --reads '/path/to/fastq/*_R{1,2}_*.fastq.gz' \
  --snp_info /path/to/info_SNPs.csv \
  --identity_amplicon AMPLxxxxxxx \
  --outdir /path/to/results \
  --threads 32 \
  -work-dir path/to/work_folder \
```

Note: Test datasets are not included in the repository due to file size.

## Example bash script for HPC

The exact submission script depends on your HPC scheduler (e.g. LSF, Slurm). The pipeline uses Nextflow profiles so that software management and cluster execution can be configured separately.

The conda profile enables the Conda environments required by individual pipeline processes, while the `lsf` profile configures Nextflow to submit each process as a separate LSF job with its own CPU and memory requirements.

For an LSF system, both profiles can be enabled with:

`-profile conda,lsf`

```
#!/bin/bash
#BSUB -J scfinded
#BSUB -o scfinded.%J.out
#BSUB -e scfinded.%J.err
#BSUB -q normal
#BSUB -n 32                 # ask for 32 CPUs
#BSUB -W 0X:00              # maximum runtime: X hours
#BSUB -G xxxxx
#BSUB -R "select[mem>128000] rusage[mem=8000]"
#BSUB -M 8000

# stop the script if a command fails, if a variable is missing, or if a piped command fails.
set -euo pipefail

module load nextflow
module load conda

cd /path/to/scFindEd-nf

nextflow run main.nf -profile conda,lsf \
  --design /path/to/design_summary.csv \
  --reads '/path/to/fastq/*_R{1,2}_*.fastq.gz' \
  --snp_info /path/to/info_SNPs.csv \
  --identity_amplicon AMPLxxxxxxx \
  --outdir /path/to/results \
  --threads 32 \
  -work-dir path/to/work_folder \

```

When the `lsf` profile is used, the LSF job shown above mainly acts as the Nextflow controller. Individual pipeline processes are submitted as separate LSF jobs, and their CPU and memory requirements are defined in `nextflow.config`. Therefore, the controller job itself does not need to request the combined resources required by the full workflow.

Cluster-specific settings such as queue names, project or group identifiers, runtime limits, and module names should be adjusted for the local HPC environment.

For execution without an LSF scheduler, the conda profile can be used on its own:

```
nextflow run main.nf -profile conda \
  [pipeline options]
```

The resolved configuration can be inspected before running with:

`nextflow config -profile conda,lsf`

This is useful for checking which executor, CPU and memory settings will be applied.

## Parameters
| Parameter   | Description                                   |
|-------------|-----------------------------------------------|
| `--design`  | Amplicon design CSV file.                      |
| `--reads`   | Paired-end FASTQ file pattern. Input FASTQ files are expected to contain cell barcodes in the read names.         |
| `--snp_info`| CSV file describing the SNPs to analyse, including HDR template sequences.         |
| `--identity_amplicon`| AmpID corresponding to the identity barcode amplicon eg: `AMPL904013,AMPL904569`.       |
| `--outdir`  | Directory where results will be written       |
| `--threads` | Number of CPU threads allocated to multithreaded processes (`default: 10`).         |
| `--compress`| Compress TSV outputs using gzip (`true` by default).        |
| `--min_mapq`| Minimum mapping quality retained after alignment (`default: 20`). Set to `0` to disable MAPQ filtering.         |


## Required inputs

1. Design Summary (`--design`)
2. SNP information (`--snp_info`)
3. FASTQ files (`--reads`)

Input FASTQ read names are expected to contain an identifier encoded using the `RG:Z:` tag, following the format: `@<read_name> RG:Z:<barcode>`. For example: `@M02088:15:000000000-LLRMH:1:1102:15466:1334 RG:Z:GGAGATGAGAAGAGAACT-1`. 

 In the current single-cell workflow, the `RG:Z:` value is treated as the cell barcode. The same tagging convention could also be used in future for sample- or experiment-level identifiers.

There is no fixed barcode length enforced by the current code. The parser simply takes everything after `RG:Z:` to the end of the FASTQ header. 

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

The `--outdir` directory contains the files required for downstream analysis with the `scFindEd` R package. Intermediate files generated during workflow execution (for example, BAM files) are not copied to the results directory by default.

These intermediate files remain in the Nextflow `work/` directory, where they are available quality control or custom downstream analyses if required. As the `work/` directory can occupy substantial disk space, it is good practice to remove it (set up in the flag `-work-dir path/to/work_folder`) once the workflow has completed successfully and the final results have been verified. Users wishing to retain intermediate files should archive or copy them before cleaning up the work/ directory.

Example structure:

```
results/
│ 
├── design/
│   └── 6915_design_summary_joined.csv
|
├── fastqc/
│   ├── sample_R1_fastqc.html
│   ├── sample_R1_fastqc.zip
│   ├── sample_R2_fastqc.html
│   └── sample_R2_fastqc.zip
│  
│ 
├── multiqc/
│   ├── multiqc_report.html
│   └── multiqc_data/
│
├── reference_amplicon/
│
├── fragment_mutations/
│   └── amplID_fragment_mutations.tsv.gz
│   └── amplID_fragment_mutations.tsv.gz
│   └── amplID_fragment_mutations.tsv.gz
│   └── amplID_fragment_mutations.tsv.gz
│ 
├── identity_amplicon/
│   └── sample_name_amplID_identity_reads.tsv.gz
│
└── qc/
     └── sample_flagstat.txt
     └── sample_idxstats.txt 
```

### Alignment filtering

After alignment with `bwa mem`, alignments are filtered before downstream processing to retain only primary mapped reads. The following alignment types are excluded:

* Unmapped reads (FLAG 4)
* Secondary alignments (FLAG 256)
* Supplementary alignments (FLAG 2048)

Filtering is performed using:

`samtools view -b -F 2308`

A default MAPQ filter of 20 is applied after alignment. This can be adjusted with `--min_mapq`, or disabled by setting `--min_mapq 0`. Disabling the MAPQ filter may be useful for highly similar amplicons where reads may be correctly aligned but assigned lower mapping quality.

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
   The SAM `FLAG` field is a single integer that stores multiple alignment properties for each read. Each property has a specific value, and the final `FLAG` is the sum of all properties that apply to that read.

   [Decode](https://broadinstitute.github.io/picard/explain-flags.html) this value using bitwise operations ([see page 7](https://samtools.github.io/hts-specs/SAMv1.pdf)) to identify which properties are present.

   The following flags are retained:
   - paired reads (`FLAG & 1`)
   - properly paired reads (`FLAG & 2`)
   - R1 reads (`FLAG & 64`)
   - R2 reads (`FLAG & 128`)

   The following flags are removed:
   - secondary alignments (`FLAG & 256`)
   - supplementary alignments (`FLAG & 2048`)

   For example:
   - `99 = 1 + 2 + 32 + 64`, meaning paired, properly paired, mate on reverse strand, and R1
   - `147 = 1 + 2 + 16 + 128`, meaning paired, properly paired, read on reverse strand, and R2

   After decoding, reads are separated into R1 and R2 groups and paired using:
   - `read_name`
   - `cell_barcode`
   - `amp_id`

   Incomplete pairs (`R1_only` and `R2_only`) could be discarded if both R1 and R2 are required to be present after pairing.

3. **Filter rows by NM values**  
   Separate rows where `r1_nm > 0` or `r2_nm > 0`, because rows with `NM = 0` for both reads have no mismatches to extract.

4. **Parse CIGAR strings**  
   The `CIGAR` string can be split into operation lengths and operation types, such as `96M`, `38S60M`, or `51M1D99M`. This can be used to track how read positions relate to amplicon-reference positions, including soft clipping, insertions and deletions.

5. **Parse MD strings**  
   The `MD` string can be used to identify mismatch positions and reference bases. For example, `73T22` means 73 matching bases, then a mismatch where the reference base is `T`, followed by 22 matching bases.

6. **Map mismatch to amplicon position**  
   The mismatch position identified from the `MD` string can be converted into an amplicon coordinate using the read alignment start position.

   For example:

   ```text
   r1_pos = 5
   r1_md  = 73T22
   ```

   means:
   - the read alignment starts at amplicon position `5`
   - the first `73` aligned bases match the reference
   - the next aligned base is a mismatch where the reference base is `T`

   Therefore, the mismatch occurs after `73` matching aligned reference positions:

   ```text
   ref_index = 74
   ```

   where `ref_index` represents the mismatch offset within the aligned reference sequence.

   The final amplicon coordinate is then calculated as:

   ```text
   amplicon_pos = r1_pos + ref_index - 1
   ```

   giving:

   ```text
   amplicon_pos = 5 + 74 - 1 = 78
   ```

   meaning the mismatch occurs at amplicon position `78`.

7. **Extract observed read base**  
   The `CIGAR` string can be used to map the mismatch reference position back to the correct position in the read sequence. The observed base can then be extracted from `r1_seq` or `r2_seq`.

8. **Create mutation call**  
   The mutation can be reported as `ref_base > read_base`, for example `T>A`.

9. **Compare with mate read**  
   For each mutation detected in R1, downstream analysis could check what base R2 has at the same amplicon position, and vice versa.

10. **Classify paired-read support**  
   Each mutation can be classified as:
   - `R1_nonoverlap_only`: mutation detected in R1, but R2 does not cover the position
   - `R2_nonoverlap_only`: mutation detected in R2, but R1 does not cover the position
   - `R1_R2_agree`: both reads cover the position and show the same mutated base
   - `R1_R2_disagree`: both reads cover the position but show different bases

11. **Propagate fragment support**  
   The `n_fragments` value can be carried forward so downstream workflows can weight each mutation by the number of identical fragments supporting that evidence pattern.

## Pipeline Reports
The pipeline does NOT automatically generates execution reports. 

## Development Notes
The pipeline is implemented in Nextflow DSL2.

Intermediate files are stored in the `work/ `directory and final results are published to the directory specified by `--outdir`.
