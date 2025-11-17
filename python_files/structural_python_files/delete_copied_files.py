
import os
import json
from pathlib import Path

def delete_copied_files(log_file="copied_files.json"):
    log_path = Path(log_file)
    if not log_path.exists():
        print(f"No log file found at {log_file}. Nothing to delete.")
        return

    with open(log_path, "r") as f:
        copied_files = json.load(f)

    deleted_count = 0
    for file_path in copied_files:
        file_path = Path(file_path)
        if file_path.exists():
            try:
                file_path.unlink()
                deleted_count += 1
            except Exception as e:
                print(f"Failed to delete {file_path}: {e}")

    # Remove log file after deletion
    log_path.unlink(missing_ok=True)

    print(f"Deleted {deleted_count} files. Log file removed.")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Delete previously copied files")
    parser.add_argument("--log", default="copied_files.json", help="Log file that tracks copied files")
    args = parser.parse_args()

    delete_copied_files(args.log)
