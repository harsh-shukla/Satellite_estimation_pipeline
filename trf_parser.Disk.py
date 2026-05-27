import subprocess
import tempfile
import os
import shutil
import sys
# from Bio import SeqIO # Removed BioPython dependency
from concurrent.futures import ProcessPoolExecutor, as_completed
import time
import math # For ceiling function
import argparse
import shlex # For parsing the options string
import gzip # For handling gzipped files

# --- Custom Parser/Writer Functions ---

def parse_fastq_gz(filename):
    """
    Generator to parse a gzipped FASTQ file and yield (header, sequence) tuples.
    Basic implementation, less robust than Bio.SeqIO for edge cases.

    Args:
        filename (str): Path to the gzipped FASTQ file.

    Yields:
        tuple: (header_line, sequence_string)
               header_line includes the starting '@'.
    """
    line_num = 0
    try:
        with gzip.open(filename, "rt", encoding='utf-8', errors='replace') as f:
            while True:
                # Read 4 lines for each record
                header = f.readline()
                if not header: # End of file
                    break
                sequence = f.readline()
                plus = f.readline()
                quality = f.readline()
                line_num += 4

                # Basic validation and cleanup
                header = header.strip()
                sequence = sequence.strip()
                plus = plus.strip()
                # quality = quality.strip() # Not needed for FASTA output

                if not header.startswith('@') or not plus.startswith('+') or not sequence:
                     # Handle potential format errors or EOF within record
                     print(f"Warning: Skipping potentially malformed FASTQ record near line {line_num} starting with: {header[:50]}...", file=sys.stderr)
                     # Attempt to find next '@' or break? For simplicity, just continue.
                     continue

                # Yield header (with '@') and sequence
                yield (header, sequence)
    except FileNotFoundError:
         print(f"Error: Input file not found during parsing: {filename}", file=sys.stderr)
         raise # Re-raise for main error handling
    except Exception as e:
         print(f"Error parsing FASTQ file {filename} near line {line_num}: {e}", file=sys.stderr)
         raise # Re-raise

def parse_fasta_gz(filename):
    """
    Generator to parse a gzipped FASTA file and yield (header, sequence) tuples.
    """
    line_num = 0
    try:
        with gzip.open(filename, "rt", encoding='utf-8', errors='replace') as f:
            header = None
            sequence_parts = []
            for line in f:
                line_num += 1
                line = line.strip()
                if not line:
                    continue
                if line.startswith('>'):
                    if header is not None:
                        yield (header, "".join(sequence_parts))
                    header = line
                    sequence_parts = []
                else:
                    if header is not None:
                        sequence_parts.append(line)
                    else:
                        print(f"Warning: Skipping sequence data without header near line {line_num}", file=sys.stderr)
            if header is not None and sequence_parts:
                yield (header, "".join(sequence_parts))
    except FileNotFoundError:
         print(f"Error: Input file not found during parsing: {filename}", file=sys.stderr)
         raise
    except Exception as e:
         print(f"Error parsing FASTA file {filename} near line {line_num}: {e}", file=sys.stderr)
         raise

def write_fasta_record(file_handle, header, sequence, line_width=60):
    """Writes a single record to an open file handle in FASTA format."""
    # Write header (replace initial '@' with '>')
    file_handle.write(f">{header[1:]}\n")
    # Write sequence with wrapping
    for i in range(0, len(sequence), line_width):
        file_handle.write(sequence[i:i+line_width] + "\n")

# --- Helper Functions ---

# Worker function remains the same as it operates on temp FASTA files
def run_trf_on_chunk_redirect(fasta_chunk_path, trf_mod_executable, trf_options_list, run_temp_dir):
    """
    Worker function executed in parallel.
    Takes a path to a temporary FASTA chunk, runs TRF-mod redirecting its stdout
    to a temporary BED file within run_temp_dir, and returns the path to that BED file.
    """
    output_bed_path = None # Initialize in case of early exit
    worker_pid = os.getpid() # Get worker PID for logging
    chunk_basename = os.path.basename(fasta_chunk_path)

    try:
        # 1. Generate a unique temporary output BED filename within run_temp_dir
        prefix = f"trf_out_{os.path.splitext(chunk_basename)[0]}_"
        temp_bed_obj = tempfile.NamedTemporaryFile(mode='w', suffix=".bed", prefix=prefix, dir=run_temp_dir, delete=False)
        output_bed_path = temp_bed_obj.name
        temp_bed_obj.close()
        if os.path.exists(output_bed_path):
             os.remove(output_bed_path) # Ensure clean slate

        # 2. Construct the TRF-mod command
        cmd = [trf_mod_executable] + trf_options_list + [fasta_chunk_path]

        # 3. Run TRF-mod, redirecting stdout to the temporary file
        stderr_output = ""
        exit_code = -1
        try:
            with open(output_bed_path, 'wb') as outfile_handle:
                process = subprocess.run(
                    cmd,
                    stdout=outfile_handle,
                    stderr=subprocess.PIPE,
                    text=False,
                    check=False
                )
            exit_code = process.returncode
            stderr_output = process.stderr.decode(errors='replace') if process.stderr else ""
        except Exception as subproc_err:
             print(f"[Worker {worker_pid}] Error during subprocess execution for {chunk_basename}: {subproc_err}", file=sys.stderr)
             raise

        # 4. Check for non-zero exit code first
        if exit_code != 0:
            error_message = (
                f"TRF-mod failed on chunk (input: {chunk_basename}) "
                f"with exit code {exit_code}.\n"
                f"Command: {' '.join(cmd)}\nStderr:\n{stderr_output}"
            )
            print(f"[Worker {worker_pid}] {error_message}", file=sys.stderr) # Log error from worker
            if output_bed_path and os.path.exists(output_bed_path):
                try: os.remove(output_bed_path) # Clean up failed output
                except OSError: pass
            raise RuntimeError(error_message) # Raise error AFTER logging stderr

        # 5. Check if output file exists and has content (even if exit code was 0)
        file_exists = os.path.exists(output_bed_path)
        file_size = os.path.getsize(output_bed_path) if file_exists else 0

        if not file_exists or file_size == 0:
            if file_exists: # Remove empty file silently
                try: os.remove(output_bed_path)
                except OSError: pass
            return None # Indicate no results

        # 6. If exit code was 0 and file has content
        return output_bed_path

    except Exception as e:
        # Catch any other exception in the worker
        if output_bed_path and os.path.exists(output_bed_path):
             try: os.remove(output_bed_path) # Clean up output file on any worker error
             except OSError: pass
        print(f"[Worker {worker_pid}] Unexpected error for chunk {chunk_basename}: {e}", file=sys.stderr)
        raise


# Merging function remains the same
def merge_temp_bed_files(temp_bed_paths, output_bed_file):
    """
    Merges temporary BED files into a single final output file.
    Assumes no headers. Deletes temp files after merging their content.
    """
    print(f"\nMerging content from temporary BED files into {output_bed_file}...")
    files_merged = 0
    lines_written = 0 # Approximate count
    valid_temp_paths = [p for p in temp_bed_paths if p is not None]
    print(f"Found {len(valid_temp_paths)} temporary BED files with content to merge.")

    try:
        # Use binary mode ('wb' and 'rb') for efficiency when just copying bytes
        with open(output_bed_file, 'wb') as outfile:
            for temp_file_path in valid_temp_paths:
                if not os.path.exists(temp_file_path):
                    print(f"Warning: Temp file {temp_file_path} not found during merge. Skipping.", file=sys.stderr)
                    continue

                try:
                    # Copy file content in chunks
                    with open(temp_file_path, 'rb') as infile:
                        while True:
                            chunk = infile.read(8192) # Read in 8KB chunks
                            if not chunk:
                                break
                            outfile.write(chunk)
                            lines_written += chunk.count(b'\n') # Approx line count

                    files_merged += 1
                    # --- CLEANUP RE-ENABLED ---
                    try:
                        os.remove(temp_file_path)
                    except OSError as e:
                        print(f"Warning: Could not remove temp BED file {temp_file_path} after merging: {e}", file=sys.stderr)
                    # --- END CLEANUP RE-ENABLED ---

                except Exception as e:
                     print(f"Warning: Failed to process temp file {temp_file_path}: {e}. Skipping.", file=sys.stderr)

        print(f"Merging complete. Merged content from {files_merged} temp files ({lines_written} lines approx).")
        if files_merged == 0 and len(valid_temp_paths) > 0:
             print("Warning: No content was merged, though temp files were expected.", file=sys.stderr)
             if lines_written == 0: open(output_bed_file, 'w').close() # Create empty file

    except IOError as e:
        print(f"Fatal Error: Could not write to final output file {output_bed_file}: {e}", file=sys.stderr)


# --- Main Orchestration Function ---

def parallel_trf_mod(input_file_gz, output_bed, num_threads, trf_mod_executable, parsed_trf_options):
    """
    Main wrapper function using a two-pass approach: count, then process.
    Reads gzipped FASTQ or FASTA input using custom parser, creates temporary FASTA chunks using custom writer.
    Workers run TRF-mod redirecting stdout to temporary BED files. Results are merged from files.
    Includes automatic cleanup of temporary files/directory. NO BIOPYTHON.
    """
    main_start_time = time.time()
    run_temp_dir = None # Initialize variable

    # --- Input Validation ---
    # (Same as before)
    if not os.path.exists(trf_mod_executable) or not os.path.isfile(trf_mod_executable):
        print(f"Error: TRF-mod executable not found or is not a file at '{trf_mod_executable}'", file=sys.stderr)
        sys.exit(1)
    if not os.access(trf_mod_executable, os.X_OK):
         print(f"Error: TRF-mod file is not executable: '{trf_mod_executable}'", file=sys.stderr)
         sys.exit(1)
    if not os.path.exists(input_file_gz):
        print(f"Error: Input file not found at '{input_file_gz}'", file=sys.stderr)
        sys.exit(1)
    if num_threads <= 0:
        print(f"Error: Number of threads must be positive.", file=sys.stderr)
        sys.exit(1)

    print(f"--- Starting Parallel TRF-mod (FASTQ/FASTA.gz Input / File BED Output / No BioPython) ---")
    print(f"Input file: {input_file_gz}")
    print(f"Output BED: {output_bed}")
    print(f"TRF-mod Path: {trf_mod_executable}")
    print(f"Requested Workers: {num_threads}")
    print(f"TRF Options: {' '.join(parsed_trf_options)}")
    print("Note: Using custom FASTQ/FASTA parser/writer (potential performance impact).")

    # Detect file type
    is_fasta = False
    lower_filename = input_file_gz.lower()
    if lower_filename.endswith(".fa.gz") or lower_filename.endswith(".fasta.gz"):
        is_fasta = True
        print("Detected FASTA format based on file extension.")
    elif lower_filename.endswith(".fq.gz") or lower_filename.endswith(".fastq.gz"):
        is_fasta = False
        print("Detected FASTQ format based on file extension.")
    else:
        print(f"Warning: Could not detect file format from extension for {input_file_gz}. Defaulting to FASTQ.", file=sys.stderr)
        is_fasta = False

    parser_func = parse_fasta_gz if is_fasta else parse_fastq_gz


    # --- Pass 1: Count Sequences and BP using Custom Parser ---
    print("\nPass 1: Counting sequences and base pairs using custom parser...")
    start_count_time = time.time()
    total_sequences = 0
    total_bp = 0
    try:
        # Use the custom generator to count records and calculate base pairs
        for _, seq in parser_func(input_file_gz):
            total_sequences += 1
            total_bp += len(seq)
    except Exception as e:
        # Error message already printed by parser, just exit
        sys.exit(1)
    count_time = time.time() - start_count_time
    print(f"Found {total_sequences} sequences and {total_bp} total base pairs in {count_time:.2f} seconds.")
    if total_sequences == 0:
        print("Error: No sequences found in input file.", file=sys.stderr)
        open(output_bed, 'w').close(); sys.exit(0)

    # Write total_bp to file
    total_bp_file = output_bed + ".total_bp"
    try:
        with open(total_bp_file, "w") as f_bp:
            f_bp.write(str(total_bp) + "\n")
        print(f"Wrote total base pairs to {total_bp_file}")
    except IOError as e:
        print(f"Warning: Could not write total_bp file: {e}", file=sys.stderr)

    # --- Calculate Chunk Size ---
    # (Same as before)
    actual_num_workers = min(num_threads, total_sequences)
    if actual_num_workers < num_threads:
        print(f"Warning: Fewer sequences ({total_sequences}) than requested threads ({num_threads}). Using {actual_num_workers} workers.")
    sequences_per_chunk = math.ceil(total_sequences / actual_num_workers)
    print(f"Calculated chunk size: {sequences_per_chunk} sequences per chunk for {actual_num_workers} target workers.")

    # --- Setup Temporary Directory and File Tracking ---
    temp_fasta_chunk_paths = []
    temp_bed_result_paths = [] # Store paths to temp BED files
    futures = []
    chunk_records = [] # Now stores (header, sequence) tuples
    sequences_in_current_chunk = 0
    chunks_submitted = 0
    start_processing_time = time.time()

    try:
        # Create the dedicated temporary directory
        run_temp_dir = tempfile.mkdtemp(prefix="trf_parallel_run_", dir='.')
        print(f"Using temporary directory for FASTA chunks and BED results: {os.path.abspath(run_temp_dir)}")

        # --- Pass 2: Process Chunks using Custom Parser and FASTA Writer ---
        print(f"\nPass 2: Processing input and submitting chunks to {actual_num_workers} workers...")
        # Use the custom parser generator
        record_iterator = parser_func(input_file_gz)
        with ProcessPoolExecutor(max_workers=actual_num_workers) as executor:
            for i, (header, sequence) in enumerate(record_iterator):
                # Collect header, sequence tuples
                chunk_records.append((header, sequence))
                sequences_in_current_chunk += 1

                is_last_sequence = (i == total_sequences - 1)
                # Check if chunk is full OR it's the very last sequence being processed
                if sequences_in_current_chunk >= sequences_per_chunk or is_last_sequence:
                    current_chunk_num = chunks_submitted + 1
                    temp_fasta_path = None

                    try:
                        # Create temporary FASTA file using custom writing logic
                        with tempfile.NamedTemporaryFile(mode='w', suffix=".fasta", prefix=f"chunk_{current_chunk_num}_", delete=False, dir=run_temp_dir) as temp_fasta:
                            temp_fasta_path = temp_fasta.name
                            # Write collected records to this file handle in FASTA format
                            for h, s in chunk_records:
                                write_fasta_record(temp_fasta, h, s) # Use custom writer
                            temp_fasta_chunk_paths.append(temp_fasta_path)

                        # Submit the job - worker function takes the temp FASTA path
                        future = executor.submit(run_trf_on_chunk_redirect, temp_fasta_path, trf_mod_executable, parsed_trf_options, run_temp_dir)
                        futures.append(future)
                        chunks_submitted += 1

                    except Exception as e_submit:
                         print(f"\n--- ERROR during preparation or submission of chunk {current_chunk_num} ---", file=sys.stderr)
                         if temp_fasta_path:
                             print(f"Attempted FASTA path: {temp_fasta_path}", file=sys.stderr)
                         print(f"Error: {e_submit}", file=sys.stderr)
                         print("Attempting to continue with other chunks...", file=sys.stderr)

                    finally:
                         # Reset for next chunk regardless of success/failure
                         chunk_records = []
                         sequences_in_current_chunk = 0

            print(f"\nFinished submitting {chunks_submitted} chunks to workers.")
            processing_submission_time = time.time() - start_processing_time
            print(f"Chunking and submission took: {processing_submission_time:.2f} seconds.")

            # --- Process Results ---
            print("\nWaiting for TRF-mod processes to complete...")
            processed_chunks = 0
            start_completion_time = time.time()
            # Collect results (paths to temp BED files or None)
            for future in as_completed(futures):
                 processed_chunks += 1
                 try:
                    # Get the result (path string or None)
                    bed_file_path_result = future.result()
                    temp_bed_result_paths.append(bed_file_path_result) # Append path or None
                 except Exception as e:
                    print(f"\n--- ERROR retrieving result for chunk (approx {processed_chunks}/{chunks_submitted}) ---", file=sys.stderr)
                    print(f"Future raised: {type(e).__name__}. See worker logs above for details.", file=sys.stderr)
                    print(f"--- END ERROR ---", file=sys.stderr)
                    temp_bed_result_paths.append(None) # Mark as failed

                 if processed_chunks % 50 == 0 or processed_chunks == chunks_submitted:
                     elapsed = time.time() - start_completion_time
                     print(f"  Processed {processed_chunks}/{chunks_submitted} chunks... ({elapsed:.1f}s elapsed)")

            completion_time = time.time() - start_completion_time
            print(f"\nFinished processing all chunks.")
            print(f"TRF-mod execution and result collection took: {completion_time:.2f} seconds.")

            # --- Combine Results from temporary BED files ---
            merge_temp_bed_files(temp_bed_result_paths, output_bed) # Use file merging function

    # Handle potential errors during the generator iteration itself
    except Exception as e:
        print(f"An unexpected error occurred during Pass 2 FASTQ processing or job submission: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1) # Ensure finally block runs
    finally:
        # --- Cleanup ---
        # Remove the entire temporary directory
        if run_temp_dir and os.path.exists(run_temp_dir):
            print(f"\nCleaning up temporary directory: {os.path.abspath(run_temp_dir)}")
            try:
                shutil.rmtree(run_temp_dir)
                print(f"Successfully removed temporary directory.")
            except Exception as cleanup_error:
                print(f"Warning: Failed to completely remove temp directory {os.path.abspath(run_temp_dir)}. Please remove it manually if needed. Error: {cleanup_error}", file=sys.stderr)
        elif temp_fasta_chunk_paths: # Fallback
             print("\nCleaning up individual temporary FASTA chunk files (fallback)...")
             for fpath in temp_fasta_chunk_paths:
                 if os.path.exists(fpath):
                     try: os.remove(fpath)
                     except OSError as e: print(f"Warning: Could not remove temp FASTA chunk {os.path.abspath(fpath)}: {e}", file=sys.stderr)


    main_end_time = time.time()
    print(f"\n--- Total execution time: {main_end_time - main_start_time:.2f} seconds ---")


# --- Command Line Argument Parsing ---
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Parallel wrapper for TRF-mod (File BED Output Mode), reading FASTQ.gz/FASTA.gz input (No BioPython), splitting based on thread count.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    # Update help text for input file
    parser.add_argument("-i", "--input", required=True, help="Input gzipped FASTQ or FASTA file path (.fastq.gz, .fq.gz, .fasta.gz, .fa.gz).")
    parser.add_argument("-o", "--output", required=True, help="Output BED file path.")
    parser.add_argument("-e", "--executable", required=True, help="Path to the TRF-mod executable.")
    parser.add_argument(
        "-t", "--threads",
        required=True,
        type=int,
        help="Number of parallel threads to use. Input FASTQ will be split into this many chunks."
    )
    parser.add_argument(
        "--trf_options",
        required=True,
        help="String containing all options for TRF-mod (e.g., \"-p 13 -s 1000\"). Quote if options contain spaces."
    )
    args = parser.parse_args()
    if args.threads <= 0:
         parser.error("Number of threads must be a positive integer.")
    try:
        parsed_trf_options = shlex.split(args.trf_options)
        print(f"Parsed TRF options: {parsed_trf_options}")
    except ValueError as e:
        parser.error(f"Error parsing --trf_options string '{args.trf_options}': {e}")

    input_file_gz = args.input
    parallel_trf_mod(
        input_file_gz=input_file_gz,
        output_bed=args.output,
        num_threads=args.threads,
        trf_mod_executable=args.executable,
        parsed_trf_options=parsed_trf_options
    )
