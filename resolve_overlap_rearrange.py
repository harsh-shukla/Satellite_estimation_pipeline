#!/usr/bin/env python3
import argparse
import sys
import os
from collections import defaultdict

def normalize_motif(pattern):
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

def resolve_overlaps(records):
    """
    Resolves overlaps within a list of records for a single contig,
    keeping the record with the largest span in case of overlap.

    Args:
        records (list): A list of tuples, where each tuple is
                        (start, end, span, original_line_string).

    Returns:
        tuple: A tuple containing two lists:
               - merged_records (list): Tuples of selected records (start, end, span, original_line).
               - discarded_lines (list): Original line strings of discarded records due to overlap.
    """
    if not records:
        return [], [] # Return two empty lists

    # Sort records primarily by start coordinate, then secondarily by end coordinate
    sorted_records = sorted(records, key=lambda x: (x[0], x[1]))

    merged_records = []
    discarded_lines = [] # New list to collect discarded lines
    if not sorted_records:
        return merged_records, discarded_lines

    # Initialize with the first record as the current best candidate
    current_best = sorted_records[0]

    for i in range(1, len(sorted_records)):
        next_record = sorted_records[i]
        # Unpack tuples for clarity
        current_start, current_end, current_span, current_line = current_best
        next_start, next_end, next_span, next_line = next_record

        # Check for overlap: next_start must be less than current_end
        if next_start < current_end:
            # Overlap exists. Compare spans and keep the longer one as current_best.
            if next_span > current_span:
                # Discard the current_best record because next_record is longer
                discarded_lines.append(current_line) # Add original line of discarded record
                current_best = next_record # Update current_best
            else:
                # Discard the next_record because current_best is longer or equal
                discarded_lines.append(next_line) # Add original line of discarded record
                # current_best remains the same
        else:
            # No overlap with the current_best record.
            # Finalize current_best for the previous group.
            merged_records.append(current_best)
            # The next_record starts a new potential group.
            current_best = next_record

    # After the loop, add the last 'current_best' record
    merged_records.append(current_best)

    return merged_records, discarded_lines # Return both lists

def main():
    parser = argparse.ArgumentParser(
        description="Filter overlapping records in a TRF-mod BED file, keeping the one with the largest span. Applies minimum length filter and normalizes motif.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument("-i", "--input", required=True, help="Input BED file path (tab-separated).")
    parser.add_argument("-o", "--output", required=True, help="Output filtered BED file path.")
    parser.add_argument(
        "--min-alignment-length",
        type=int,
        default=0, # Default to 0 (no filtering) if not provided
        help="Minimum span (end - start) required for a record to be kept in the final output."
        )
    args = parser.parse_args()

    # Define the discard file name
    discard_filename = "Discarded_overlaps.txt"

    if not os.path.exists(args.input):
        print(f"Error: Input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    if args.min_alignment_length < 0:
         print(f"Error: --min-alignment-length cannot be negative.", file=sys.stderr)
         sys.exit(1)

    # Use defaultdict for easier handling of contigs
    contig_data = defaultdict(list)
    header_lines = []
    print(f"Reading input BED file: {args.input}...")
    line_num = 0
    records_read = 0
    malformed_skipped = 0
    negative_span_skipped = 0

    try:
        with open(args.input, 'r') as infile:
            for line in infile:
                line_num += 1
                line_content = line.strip() # Process stripped content
                original_line = line # Keep original line with newline for output/discard

                if not line_content:
                    continue # Skip empty lines

                # Preserve header lines
                if line_content.startswith(('#', 'track', 'browser')):
                    header_lines.append(original_line) # Keep original line ending
                    continue

                fields = line_content.split('\t')
                # Expecting ctg, start, end, ..., pattern (10 fields total, index 0-9)
                # *** Corrected length check ***
                if len(fields) < 10:
                    print(f"Warning: Skipping malformed line {line_num} (less than 10 fields): '{line_content}'", file=sys.stderr)
                    malformed_skipped += 1
                    continue

                try:
                    ctg = fields[0]
                    start = int(fields[1])
                    end = int(fields[2])
                    span = end - start

                    if start < 0 or end < 0:
                         print(f"Warning: Skipping record with negative coordinate on line {line_num}: '{line_content}'", file=sys.stderr)
                         malformed_skipped += 1
                         continue
                    if span < 0:
                         print(f"Warning: Skipping record with negative span (end < start) on line {line_num}: '{line_content}'", file=sys.stderr)
                         negative_span_skipped += 1
                         continue

                    # Store relevant info: start, end, span, and the original full line
                    record_tuple = (start, end, span, original_line)

                    contig_data[ctg].append(record_tuple)
                    records_read += 1

                except ValueError:
                    print(f"Warning: Skipping line {line_num} due to non-integer start/end fields: '{line_content}'", file=sys.stderr)
                    malformed_skipped += 1
                    continue
                except IndexError: # Should be caught by len check now
                     print(f"Warning: Skipping severely malformed line {line_num}: '{line_content}'", file=sys.stderr)
                     malformed_skipped += 1
                     continue

    except Exception as e:
         print(f"Error reading input file {args.input}: {e}", file=sys.stderr)
         sys.exit(1)

    print(f"Read {records_read} valid data records across {len(contig_data)} contigs.")
    if malformed_skipped > 0:
        print(f"Skipped {malformed_skipped} malformed lines.")
    if negative_span_skipped > 0:
         print(f"Skipped {negative_span_skipped} records with negative span.")

    print("Resolving overlaps by keeping the longest span...")
    print(f"Applying minimum alignment length filter: span >= {args.min_alignment_length}")
    print(f"Discarded overlap records will be written to: {discard_filename}")

    final_filtered_count = 0
    all_discarded_overlap_lines = [] # Collect all discarded lines here

    try:
        # Open the main output file
        with open(args.output, 'w') as outfile:
            # Write preserved header lines first
            for header in header_lines:
                outfile.write(header)

            # Process each contig
            sorted_contigs = sorted(contig_data.keys())
            for ctg in sorted_contigs:
                records = contig_data[ctg]
                # Resolve overlaps first
                resolved_records, discarded_overlap_lines = resolve_overlaps(records)
                # Add lines discarded due to overlap to the main list
                all_discarded_overlap_lines.extend(discarded_overlap_lines)

                # Apply length filter and normalize motif for records kept after overlap resolution
                for record_tuple in resolved_records:
                    start, end, span, original_line = record_tuple

                    # Apply minimum length filter
                    if span >= args.min_alignment_length:
                        try:
                            # Parse the original line to modify the pattern
                            fields = original_line.strip().split('\t')
                            # *** Corrected length check and index ***
                            if len(fields) >= 10: # Need 10 fields to access index 9
                                pattern = fields[9] # Get pattern from index 9
                                normalized_pattern = normalize_motif(pattern)
                                # Replace the old pattern (index 9) with the normalized one
                                fields[9] = normalized_pattern
                                # Reconstruct the line
                                modified_line = "\t".join(fields) + "\n"
                                outfile.write(modified_line)
                                final_filtered_count += 1
                            else:
                                # Handle case where line somehow passed initial check but fails here
                                print(f"Warning: Could not process record for motif normalization (unexpected field count): {original_line.strip()}", file=sys.stderr)
                                outfile.write(original_line) # Write original line as fallback
                                final_filtered_count += 1

                        except Exception as e_norm:
                             print(f"Warning: Error normalizing motif for record, writing original: {original_line.strip()}. Error: {e_norm}", file=sys.stderr)
                             outfile.write(original_line) # Write original line as fallback
                             final_filtered_count += 1
                    # else: record is discarded due to length filter (do nothing further with it)


        # After processing all contigs, write the discarded overlap records to their file
        print(f"\nWriting {len(all_discarded_overlap_lines)} discarded overlap records to {discard_filename}...")
        with open(discard_filename, 'w') as discard_file:
            discard_file.write("# Records discarded during overlap resolution (longest span kept)\n")
            for line in all_discarded_overlap_lines:
                 # Write original line prefixed, stripping its newline and adding one
                 discard_file.write(f"# Discarded:\t{line.strip()}\n")

    except IOError as e:
         print(f"Error writing to output file {args.output} or {discard_filename}: {e}", file=sys.stderr)
         sys.exit(1)

    print(f"\nFinished.")
    print(f"Wrote {final_filtered_count} selected and filtered records to {args.output}.")
    print(f"Wrote {len(all_discarded_overlap_lines)} discarded overlap records to {discard_filename}.")

if __name__ == "__main__":
    main()
