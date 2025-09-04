import os
import sys

def get_sample_paths(sample_list_file):
    valid_samples = []
    invalid_samples = []

    with open(sample_list_file) as sample_list:
        next(sample_list)  # Skip header
        for line_number, line in enumerate(sample_list, start=2):
            line = line.strip().split('\t')
            if len(line) >= 3:
                sample = line[0].strip()
                path = line[2].strip()

                if os.path.isdir(path):
                    print(f"[PASS] Sample '{sample}' ({sample_list_file} line {line_number}): Directory exists -> {path}", file=sys.stderr)
                    valid_samples.append(path)
                else:
                    print(f"[FAIL] Sample '{sample}' ({sample_list_file} line {line_number}): Directory does NOT exist -> {path}", file=sys.stderr)
                    invalid_samples.append((sample, line_number, path))
            else:
                print(f"[ERROR] {sample_list_file} line {line_number}: Incomplete or malformed entry.", file=sys.stderr)

    if invalid_samples:
        print("\n=== INVALID PATHS DETECTED ===", file=sys.stderr)
        for sample, line_number, path in invalid_samples:
            print(f"Sample '{sample}' ({sample_list_file} line {line_number}): Invalid path -> {path}", file=sys.stderr)
        sys.exit(f"\nExiting: {len(invalid_samples)} invalid path(s) found.")

    # Print comma-sep paths on one line (to stdout)
    print(",".join(valid_samples))

    return valid_samples

if __name__ == "__main__":
    get_sample_paths('samples.sample_list')
