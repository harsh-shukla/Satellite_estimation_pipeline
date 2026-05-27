#!/usr/bin/env python3
import argparse
import sys
import os
from collections import defaultdict

# --- Helper Functions ---

def reverse_complement(dna):
    """Computes the reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
                  'N': 'N', # Handle Ns if present
                  'a': 't', 't': 'a', 'c': 'g', 'g': 'c', 'n': 'n'}
    try:
        # Reverse the sequence and complement each base
        return "".join(complement[base] for base in reversed(dna))
    except KeyError as e:
        print(f"Warning: Invalid character '{e}' found in sequence '{dna}'. Cannot compute reverse complement accurately.", file=sys.stderr)
        return dna # Return original on error

def get_smallest_rotation(pattern):
    """Finds the lexicographically smallest rotation of a pattern string."""
    n = len(pattern)
    if n == 0:
        return ""
    # Concatenate the pattern with itself to easily access all rotations
    double_pattern = pattern + pattern
    # Find the starting index of the lexicographically smallest rotation
    min_rotation_start = 0
    for i in range(1, n):
        # Compare the rotation starting at index i with the current minimum
        if double_pattern[i : i + n] < double_pattern[min_rotation_start : min_rotation_start + n]:
            min_rotation_start = i
    # Return the smallest rotation found
    return double_pattern[min_rotation_start : min_rotation_start + n]

def get_canonical_motif(motif):
    """
    Determines the true canonical representation of a motif by considering
    both rotations and reverse complements, returning the lexicographically
    smallest representation overall.

    Args:
        motif (str): The motif string, assumed to be the smallest rotation
                     of its strand (output from previous script).

    Returns:
        str: The true canonical motif representation.
    """
    if not motif:
        return ""

    # 1. The input motif is assumed to be its strand's smallest rotation.
    normalized_motif = motif

    # 2. Calculate its reverse complement.
    rc_motif = reverse_complement(motif)

    # 3. Find the smallest rotation *of the reverse complement*.
    normalized_rc_motif = get_smallest_rotation(rc_motif)

    # 4. The true canonical form is the smaller of the two normalized strings.
    return min(normalized_motif, normalized_rc_motif)

# --- Main Function ---

def main():
    parser = argparse.ArgumentParser(
        description="Summarize satellite repeat content from a filtered TRF BED file, using true canonical motifs.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input filtered, motif-normalized BED file path (output from resolve_bed_overlaps.py)."
        )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output summary file path (tab-separated)."
        )
    parser.add_argument(
        "--total-bp-sequenced",
        type=str,
        required=True,
        help="Total number of base pairs in the original sequencing dataset, or a file containing this number."
        )
    args = parser.parse_args()

    # --- Input Validation ---
    if not os.path.exists(args.input):
        print(f"Error: Input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)
    total_bp_seq_val = -1
    if os.path.exists(args.total_bp_sequenced) and os.path.isfile(args.total_bp_sequenced):
        try:
            with open(args.total_bp_sequenced, 'r') as f:
                content = f.read().strip()
                total_bp_seq_val = int(content)
        except Exception as e:
            print(f"Error reading total base pairs from file {args.total_bp_sequenced}: {e}", file=sys.stderr)
            sys.exit(1)
    else:
        try:
            total_bp_seq_val = int(args.total_bp_sequenced)
        except ValueError:
            print(f"Error: --total-bp-sequenced must be a positive integer or a valid file path.", file=sys.stderr)
            sys.exit(1)

    if total_bp_seq_val <= 0:
         print(f"Error: --total-bp-sequenced must be a positive integer.", file=sys.stderr)
         sys.exit(1)

    # --- Data Aggregation Structure ---
    aggregated_data = defaultdict(lambda: {
        "motif_bp": 0,
        "revcomp_bp": 0,
        "motif_contigs": set(),
        "revcomp_contigs": set()
    })

    print(f"Reading filtered BED file: {args.input}...")
    line_num = 0
    records_processed = 0
    malformed_skipped = 0

    try:
        with open(args.input, 'r') as infile:
            for line in infile:
                line_num += 1
                line_content = line.strip()
                if not line_content or line_content.startswith(('#', 'track', 'browser')):
                    continue # Skip empty lines and headers

                fields = line_content.split('\t')
                # Expecting ctg, start, end, ..., pattern (10 fields total, index 0-9)
                if len(fields) < 10:
                    print(f"Warning: Skipping malformed line {line_num} (less than 10 fields): '{line_content}'", file=sys.stderr)
                    malformed_skipped += 1
                    continue

                try:
                    ctg = fields[0]
                    start = int(fields[1])
                    end = int(fields[2])
                    # Pattern read from file - assumed to be smallest rotation of its strand
                    pattern_from_file = fields[9]
                    span = end - start

                    if span < 0:
                         print(f"Warning: Skipping record with negative span on line {line_num}: '{line_content}'", file=sys.stderr)
                         malformed_skipped += 1
                         continue
                    if not pattern_from_file:
                         print(f"Warning: Skipping record with empty pattern on line {line_num}: '{line_content}'", file=sys.stderr)
                         malformed_skipped += 1
                         continue

                    # Determine the *true* canonical representation using the corrected function
                    canonical_motif = get_canonical_motif(pattern_from_file)
                    data = aggregated_data[canonical_motif]

                    # Check if the pattern from the file matches the final canonical form.
                    # If it doesn't match, it must belong to the reverse complement group relative
                    # to the final canonical form.
                    if pattern_from_file == canonical_motif:
                        # This record's pattern *is* the final canonical form
                        data["motif_bp"] += span
                        data["motif_contigs"].add(ctg)
                    else:
                        # This record's pattern corresponds to the reverse complement
                        # of the final canonical form.
                        data["revcomp_bp"] += span
                        data["revcomp_contigs"].add(ctg)

                    records_processed += 1

                except ValueError:
                    print(f"Warning: Skipping line {line_num} due to non-integer start/end fields: '{line_content}'", file=sys.stderr)
                    malformed_skipped += 1
                    continue
                except IndexError:
                     print(f"Warning: Skipping severely malformed line {line_num}: '{line_content}'", file=sys.stderr)
                     malformed_skipped += 1
                     continue

    except Exception as e:
         print(f"Error reading input file {args.input}: {e}", file=sys.stderr)
         sys.exit(1)

    print(f"Processed {records_processed} data records.")
    if malformed_skipped > 0:
        print(f"Skipped {malformed_skipped} malformed lines.")
    print(f"Found {len(aggregated_data)} unique canonical motifs (considering rotations and strand).")

    # --- Output Generation ---
    print(f"Calculating statistics and writing to {args.output}...")
    try:
        with open(args.output, 'w') as outfile:
            # Write header
            header = "\t".join([
                "CanonicalMotif", "MotifBpCov", "MotifRatio",
                "RevCompMotif", "RevCompBpCov", "RevCompRatio",
                "TotalBpCombined", "TotalRatioCombined",
                "MotifContigs", "RevCompContigs"
            ])
            outfile.write(header + "\n")

            # Process aggregated data, sorting by canonical motif
            for motif in sorted(aggregated_data.keys()): # motif is now the true canonical form
                data = aggregated_data[motif]

                # Calculate values
                motif_bp = data["motif_bp"] # BP for the canonical strand
                revcomp_bp = data["revcomp_bp"] # BP for the reverse complement strand
                total_bp_combined = motif_bp + revcomp_bp

                # Calculate ratios
                total_bp_seq = float(total_bp_seq_val)
                motif_ratio = motif_bp / total_bp_seq if total_bp_seq > 0 else 0.0
                revcomp_ratio = revcomp_bp / total_bp_seq if total_bp_seq > 0 else 0.0
                total_ratio_combined = total_bp_combined / total_bp_seq if total_bp_seq > 0 else 0.0

                # Get the actual reverse complement string of the canonical motif
                # Note: This might not be the smallest rotation of the RC strand,
                # but it's the direct RC of the canonical motif listed in column 1.
                revcomp_motif_str = reverse_complement(motif)

                # Format contig lists
                motif_contigs_str = ",".join(sorted(list(data["motif_contigs"])))
                revcomp_contigs_str = ",".join(sorted(list(data["revcomp_contigs"])))

                # Format output line
                output_fields = [
                    motif, # The true canonical motif
                    motif_bp,
                    f"{motif_ratio:.6f}",
                    revcomp_motif_str, # The direct RC of the canonical motif
                    revcomp_bp,
                    f"{revcomp_ratio:.6f}",
                    total_bp_combined,
                    f"{total_ratio_combined:.6f}",
                    motif_contigs_str if motif_contigs_str else "NA",
                    revcomp_contigs_str if revcomp_contigs_str else "NA"
                ]
                outfile.write("\t".join(map(str, output_fields)) + "\n")

    except IOError as e:
         print(f"Error writing to output file {args.output}: {e}", file=sys.stderr)
         sys.exit(1)
    except Exception as e:
         print(f"An unexpected error occurred during output generation: {e}", file=sys.stderr)
         import traceback
         traceback.print_exc()
         sys.exit(1)

    print(f"\nFinished summarizing satellite content.")
    print(f"Summary written to {args.output}")

if __name__ == "__main__":
    main()
