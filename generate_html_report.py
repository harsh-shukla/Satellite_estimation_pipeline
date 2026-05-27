#!/usr/bin/env python3
import csv
import sys
import os
import argparse

# Increase CSV field size limit for large contig strings
maxInt = sys.maxsize
while True:
    try:
        csv.field_size_limit(maxInt)
        break
    except OverflowError:
        maxInt = int(maxInt/10)

def main():
    parser = argparse.ArgumentParser(description="Generate interactive HTML report from satellite summary TSV.")
    parser.add_argument("-i", "--input", required=True, help="Input satellite_summary.tsv file.")
    parser.add_argument("-o", "--output", default="satellite_report.html", help="Output HTML file path.")
    parser.add_argument("--min-ratio", type=float, default=0.00001, help="Minimum TotalRatioCombined to display.")
    args = parser.parse_args()

    if not os.path.exists(args.input):
        print(f"Error: Input file {args.input} not found.", file=sys.stderr)
        sys.exit(1)

    data = []
    header = []
    try:
        with open(args.input, 'r', encoding='utf-8') as f:
            reader = csv.reader(f, delimiter='\t')
            header = next(reader)
            
            # Find index of TotalRatioCombined
            try:
                ratio_idx = header.index("TotalRatioCombined")
            except ValueError:
                print("Error: 'TotalRatioCombined' column not found in TSV.", file=sys.stderr)
                sys.exit(1)
                
            # Filter out MotifContigs and RevCompContigs
            exclude_cols = ["MotifContigs", "RevCompContigs"]
            keep_indices = [i for i, col in enumerate(header) if col not in exclude_cols]
            
            header = [header[i] for i in keep_indices]
            new_ratio_idx = header.index("TotalRatioCombined")
                
            for row in reader:
                if len(row) <= ratio_idx:
                    continue
                try:
                    ratio = float(row[ratio_idx])
                except ValueError:
                    continue
                
                if ratio >= args.min_ratio:
                    filtered_row = [row[i] for i in keep_indices if i < len(row)]
                    data.append(filtered_row)
                    
    except Exception as e:
        print(f"Error reading {args.input}: {e}", file=sys.stderr)
        sys.exit(1)
        
    # Sort data descending by TotalRatioCombined
    data.sort(key=lambda x: float(x[new_ratio_idx]), reverse=True)

    # Generate HTML
    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Satellite Summary Report</title>
    <link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/1.13.6/css/jquery.dataTables.min.css">
    <style>
        body {{ font-family: sans-serif; margin: 20px; background-color: #fafafa; color: #333; }}
        h1 {{ text-align: center; color: #2c3e50; margin-bottom: 5px; }}
        .subtitle {{ text-align: center; color: #7f8c8d; margin-top: 0; margin-bottom: 30px; }}
        .container {{ width: 95%; margin: 0 auto; background: white; padding: 20px; border-radius: 8px; box-shadow: 0 4px 6px rgba(0,0,0,0.1); }}
        table {{ width: 100%; }}
        td {{ word-break: break-all; font-size: 14px; }}
        th {{ background-color: #34495e; color: white; padding: 10px; font-size: 15px; }}
    </style>
</head>
<body>
    <div class="container">
        <h1>Satellite Estimation Summary</h1>
        <p class="subtitle">Filtered for TotalRatioCombined &ge; {args.min_ratio}</p>
        <table id="satelliteTable" class="display cell-border hover row-border">
            <thead>
                <tr>
                    {"".join(f"<th>{col}</th>" for col in header)}
                </tr>
            </thead>
            <tbody>
"""
    
    for row in data:
        html_content += "                <tr>\n"
        for item in row:
            html_content += f"                    <td>{item}</td>\n"
        html_content += "                </tr>\n"
        
    html_content += """            </tbody>
        </table>
    </div>

    <script src="https://code.jquery.com/jquery-3.7.0.min.js"></script>
    <script src="https://cdn.datatables.net/1.13.6/js/jquery.dataTables.min.js"></script>
    <script>
        $(document).ready(function() {
            $('#satelliteTable').DataTable({
                "pageLength": 50,
                "order": [] // Disable initial sorting as it's already sorted in python
            });
        });
    </script>
</body>
</html>
"""

    try:
        with open(args.output, 'w', encoding='utf-8') as f:
            f.write(html_content)
        print(f"Interactive HTML report generated successfully at: {args.output}")
    except Exception as e:
        print(f"Error writing to {args.output}: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
