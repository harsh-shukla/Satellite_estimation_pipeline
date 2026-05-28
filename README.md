# Satellite Estimation Pipeline

This repository contains a parallelized pipeline for estimating satellite and tandem repeat sequences from long-read sequencing datasets (e.g., PacBio HiFi or Oxford Nanopore). It acts as a highly efficient wrapper around `TRF-mod` (Tandem Repeats Finder), capable of utilizing multiple CPU cores and rapidly processing large `FASTQ.gz` and `FASTA.gz` files.

## Prerequisites

If running locally:
* **TRF-mod**: You must have the [TRF-mod](https://github.com/lh3/TRF-mod) executable installed on your system or accessible via your HPC.
* **Python 3**: The scripts require standard Python 3. No external libraries (like `pandas` or `Biopython`) are required. This ensures the pipeline remains lightweight and operates out-of-the-box even on restrictive clusters.

Alternatively, you can run the pipeline without installing any prerequisites by using the provided **Docker** or **Singularity** container (see the Container Usage section below).

## Pipeline Components

The pipeline consists of four Python scripts and a Bash wrapper to tie them together:

### 1. `trf_parser.Disk.py`
A parallel wrapper for TRF-mod.
- **Input**: A `FASTQ.gz` or `FASTA.gz` file. The script automatically detects the format from the extension.
- **Process**: Splits the multi-gigabyte sequence file into temporary chunks based on your requested number of threads. It calculates the total sequenced base pairs while streaming the files. Then, it spins up multiple background workers to run `TRF-mod` on each chunk simultaneously.
- **Output**: A concatenated `satellite.bed` file containing all raw repeat instances across the sequenced reads. It also drops a small `satellite.bed.total_bp` file containing the overall sequenced bases.

### 2. `resolve_overlap_rearrange.py`
A cleaner and normalizer.
- **Input**: The raw `satellite.bed` file.
- **Process**: Parses the BED format and identifies any overlapping tandem repeats on the same read (a common issue with raw TRF outputs). If an overlap occurs, it keeps the motif span that is longest and discards the other. It also applies an optional minimum alignment length filter. It also normalize the motif i.e finds the lexicographically smallest rotation of a pattern string on the sequenced strand.
- **Output**: A normalized, non-overlapping BED file (e.g., `satellite.FILTNORM.bed`) and a discard log (`Discarded_overlaps.txt`).

### 3. `summarize_trf_content.py`
A summary generator.
- **Input**: The filtered BED file and the `total_bp` metric (as a raw integer or by passing the `.total_bp` file generated in step 1).
- **Process**: Groups motifs w.r.t to their canonical representation and tallies the total base pairs occupied by each satellite across all contigs. 
- **Output**: A detailed summary file (`satellite_summary.tsv`) displaying base pair coverage and the genome-wide percentage ratio of each satellite.

### 4. `generate_html_report.py`
An HTML reporting tool.
- **Input**: The `satellite_summary.tsv` file.
- **Process**: Filters out extremely low abundance satellites (by default, `< 0.00001` or 0.001%) to declutter the final view, and sorts the remaining satellites by abundance. It generates a standalone HTML document encapsulating the data.
- **Output**: An interactive `satellite_report.html` page embedded with DataTables, enabling searching, sorting, and pagination directly in any modern web browser.

## Standalone Usage

You can run the full pipeline sequentially using the provided `run_wrapper.sh` script. The wrapper script has been updated to accept command-line arguments to easily configure the underlying python scripts.

```bash
./run_wrapper.sh -i input.fastq.gz [options]

Options:
  -i, --input            Input sequence file (FASTQ/FASTA, gzipped)
  -o, --prefix           Prefix for all output files (default: satellite)
  -t, --threads          Number of threads to use for TRF (default: 1)
  --trf-options          Options string to pass to TRF-mod (default: '-p 15 -s 100')
  --trf-bin              Path to TRF-mod executable (default: 'trf-mod' assuming in PATH)
  --min-align            Minimum alignment length for resolving overlaps (default: 50)
  --min-ratio            Minimum TotalRatioCombined to display in HTML report (default: 0.00001)
```

Example workflow:
```bash
./run_wrapper.sh -i data/ISO1_hifi_099.fastq.gz -t 32 --min-align 50 --trf-bin /path/to/trf-mod
```

## Docker / Singularity Container Usage

A Docker image containing Python, TRF-mod, and all the scripts is automatically built and hosted on the GitHub Container Registry (GHCR). This is highly recommended for reproducibility and running on HPC environments.

### Using Docker

If you have Docker installed, you can run the pipeline directly. Remember to mount your data directory using `-v` so the container can access your input files and write the outputs.

```bash
# Example using Docker
docker run --rm -v /path/to/your/data:/data \
  ghcr.io/harsh-shukla/satellite_estimation_pipeline:latest \
  -i /data/input.fastq.gz \
  -o /data/satellite_output \
  -t 32
```

### Using Singularity / Apptainer (HPC Environments)

Singularity can seamlessly pull and convert Docker images from GHCR. This is the preferred method on most computing clusters.

```bash
# Pull the image and convert it to a Singularity Image Format (.sif)
singularity build satellite_pipeline.sif docker://ghcr.io/harsh-shukla/satellite_estimation_pipeline:latest

# Run the pipeline
# Singularity automatically binds your current directory and home directory by default
singularity run satellite_pipeline.sif \
  -i my_input.fastq.gz \
  -o my_satellite_output \
  -t 32 \
  --min-align 50
```

## Outputs

After a complete, successful run, your directory will contain:
* `[prefix].bed`: Raw concatenated TRF output.
* `[prefix].bed.total_bp`: Text file containing the overall total base pairs sequenced in the input dataset.
* `[prefix].FILTNORM.bed`: Filtered BED file with overlaps resolved.
* `Discarded_overlaps.txt`: Log of records discarded during overlap resolution.
* `[prefix]_summary.tsv`: Tab-separated summary file detailing the abundance of each canonical motif.
* **`[prefix]_report.html`**: An interactive HTML visualization of the summary file for easy inspection and reporting.
