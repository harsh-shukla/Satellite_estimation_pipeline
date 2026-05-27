# Satellite Estimation Pipeline

This repository contains a parallelized pipeline for estimating satellite and tandem repeat sequences from long-read sequencing datasets (e.g., PacBio HiFi or Oxford Nanopore). It acts as a highly efficient wrapper around `TRF-mod` (Tandem Repeats Finder), capable of utilizing multiple CPU cores and rapidly processing large `FASTQ.gz` and `FASTA.gz` files.

## Prerequisites

* **TRF-mod**: You must have the [TRF-mod](https://github.com/lh3/TRF-mod) executable installed on your system or accessible via your HPC.
* **Python 3**: The scripts require standard Python 3. No external libraries (like `pandas` or `Biopython`) are required. This ensures the pipeline remains lightweight and operates out-of-the-box even on restrictive clusters.

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

## Usage Guide

You can run the full pipeline sequentially using the provided `run_wrapper.sh` script. Adjust the paths to point to your specific FASTQ/FASTA data and your `TRF-mod` binary location.

Example workflow (as seen in `run_wrapper.sh`):

```bash
# 1. Run parallel TRF-mod (Modify -t to your desired thread count)
python3 trf_parser.Disk.py -i input.fastq.gz -o satellite.bed -e /path/to/trf-mod -t 32 --trf_options "-p 15 -s 100"

# 2. Resolve overlaps (minimum alignment length of 50bp)
python3 resolve_overlap_rearrange.py -i satellite.bed -o satellite.FILTNORM.bed --min-alignment-length 50

# 3. Summarize content
# Note: The total-bp-sequenced is read dynamically from the file outputted by step 1
python3 summarize_trf_content.py -i satellite.FILTNORM.bed -o satellite_summary.tsv --total-bp-sequenced satellite.bed.total_bp

# 4. Generate Interactive HTML report
python3 generate_html_report.py -i satellite_summary.tsv -o satellite_report.html
```

## Outputs

After a complete, successful run, your directory will contain:
* `satellite.bed`: Raw concatenated TRF output.
* `satellite.bed.total_bp`: Text file containing the overall total base pairs sequenced in the input dataset.
* `satellite.FILTNORM.bed`: Filtered BED file with overlaps resolved.
* `Discarded_overlaps.txt`: Log of records discarded during overlap resolution.
* `satellite_summary.tsv`: Tab-separated summary file detailing the abundance of each canonical motif.
* **`satellite_report.html`**: An interactive HTML visualization of the summary file for easy inspection and reporting.
