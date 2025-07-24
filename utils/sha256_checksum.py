import hashlib
import os
from pathlib import Path

def compute_sha256(file_path):
    """Compute SHA-256 checksum for a file."""
    sha256 = hashlib.sha256()
    with open(file_path, "rb") as f:
        while chunk := f.read(8192):
            sha256.update(chunk)
    return sha256.hexdigest()

def generate_checksums(data_dir, output_file):
    """Generate SHA-256 checksums for all files in a directory."""
    data_dir = Path(data_dir)
    output_file = Path(output_file)
    
    with open(output_file, "w") as f:
        for file in data_dir.glob("*"):
            if file.is_file() and file.name != output_file.name:  # Exclude checksums.txt
                checksum = compute_sha256(file)
                f.write(f"{file.name} {checksum}\n")
    print(f"Checksums saved to {output_file}")

def verify_checksums(data_dir, checksum_file):
    """Verify SHA-256 checksums for all files in a directory."""
    data_dir = Path(data_dir)
    checksum_file = Path(checksum_file)
    
    with open(checksum_file, "r") as f:
        saved_checksums = {line.split()[0]: line.split()[1] for line in f}
    
    for file in data_dir.glob("*"):
        if file.is_file() and file.name != checksum_file.name:  # Exclude checksums.txt
            computed_checksum = compute_sha256(file)
            saved_checksum = saved_checksums.get(file.name)
            if not saved_checksum:
                print(f"WARNING: No checksum found for {file.name}")
            elif computed_checksum != saved_checksum:
                print(f"ERROR: Checksum mismatch for {file.name}")
            else:
                print(f"Verified: {file.name}")
    print("Checksum verification complete.")

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Generate or verify SHA-256 checksums.")
    parser.add_argument("action", choices=["generate", "verify"], help="Action to perform: generate or verify checksums.")
    parser.add_argument("data_dir", help="Directory containing the data files.")
    parser.add_argument("checksum_file", help="File to save or read checksums.")
    args = parser.parse_args()

    if args.action == "generate":
        generate_checksums(args.data_dir, args.checksum_file)
    elif args.action == "verify":
        verify_checksums(args.data_dir, args.checksum_file)