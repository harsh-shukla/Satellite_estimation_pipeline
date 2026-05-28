#!/bin/bash
set -e

# Resolve script directory so python scripts can be found relative to this wrapper
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# Default arguments
THREADS=1
TRF_OPTIONS="-p 15 -s 100"
MIN_ALIGN=50
MIN_RATIO=0.00001
PREFIX="satellite"
TRF_BIN="trf-mod"

usage() {
  echo "Satellite Estimation Pipeline Wrapper"
  echo "Usage: $0 -i <input.fastq.gz|input.fasta.gz> [options]"
  echo ""
  echo "Options:"
  echo "  -i, --input            Input sequence file (FASTQ/FASTA, gzipped)"
  echo "  -o, --prefix           Prefix for all output files (default: satellite)"
  echo "  -t, --threads          Number of threads to use for TRF (default: 1)"
  echo "  --trf-options          Options string to pass to TRF-mod (default: '-p 15 -s 100')"
  echo "  --trf-bin              Path to TRF-mod executable (default: 'trf-mod' assuming in PATH)"
  echo "  --min-align            Minimum alignment length for resolving overlaps (default: 50)"
  echo "  --min-ratio            Minimum TotalRatioCombined to display in HTML report (default: 0.00001)"
  echo "  -h, --help             Show this help message and exit"
  exit 1
}

# Parse command line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -i|--input) INPUT="$2"; shift ;;
        -o|--prefix) PREFIX="$2"; shift ;;
        -t|--threads) THREADS="$2"; shift ;;
        --trf-options) TRF_OPTIONS="$2"; shift ;;
        --trf-bin) TRF_BIN="$2"; shift ;;
        --min-align) MIN_ALIGN="$2"; shift ;;
        --min-ratio) MIN_RATIO="$2"; shift ;;
        -h|--help) usage ;;
        *) echo "Unknown parameter passed: $1"; usage ;;
    esac
    shift
done

if [ -z "$INPUT" ]; then
    echo "Error: Input file (-i) is required."
    usage
fi

echo "=========================================="
echo " Starting Satellite Estimation Pipeline"
echo "=========================================="
echo "Input: $INPUT"
echo "Prefix: $PREFIX"
echo "Threads: $THREADS"
echo "TRF Options: $TRF_OPTIONS"
echo "Min Alignment: $MIN_ALIGN"
echo "Min Ratio (HTML): $MIN_RATIO"
echo "=========================================="

echo "[1/4] Running TRF-mod parallel wrapper..."
python3 "${SCRIPT_DIR}/trf_parser.Disk.py" \
    -i "$INPUT" \
    -o "${PREFIX}.bed" \
    -e "$TRF_BIN" \
    -t "$THREADS" \
    --trf_options "$TRF_OPTIONS"

echo "[2/4] Resolving overlaps..."
python3 "${SCRIPT_DIR}/resolve_overlap_rearrange.py" \
    -i "${PREFIX}.bed" \
    -o "${PREFIX}.FILTNORM.bed" \
    --min-alignment-length "$MIN_ALIGN"

echo "[3/4] Summarizing TRF content..."
python3 "${SCRIPT_DIR}/summarize_trf_content.py" \
    -i "${PREFIX}.FILTNORM.bed" \
    -o "${PREFIX}_summary.tsv" \
    --total-bp-sequenced "${PREFIX}.bed.total_bp"

echo "[4/4] Generating HTML report..."
python3 "${SCRIPT_DIR}/generate_html_report.py" \
    -i "${PREFIX}_summary.tsv" \
    -o "${PREFIX}_report.html" \
    --min-ratio "$MIN_RATIO"

echo "=========================================="
echo " Pipeline Complete!"
echo " Outputs saved with prefix: ${PREFIX}"
echo "=========================================="
