import os
import shutil
import json
from pathlib import Path

def copy_txt_files(src_folder, dest_folder, dest_2d_folder, log_file="copied_files.json"):
    # Load existing log if available
    log_path = Path(log_file)
    if log_path.exists():
        with open(log_path, "r") as f:
            copied_files = set(json.load(f))
    else:
        copied_files = set()

    src_folder = Path(src_folder).resolve()
    dest_folder = Path(dest_folder).resolve()
    dest_2d_folder = Path(dest_2d_folder).resolve()

    dest_folder.mkdir(parents=True, exist_ok=True)
    dest_2d_folder.mkdir(parents=True, exist_ok=True)

    new_files = []

    for root, _, files in os.walk(src_folder):
        for file in files:
            if file.lower().endswith(".txt"):
                source_path = Path(root) / file

                # Choose destination base folder
                if "2D" in file:
                    target_path = dest_2d_folder / file
                else:
                    target_path = dest_folder / file

                # Skip if file already exists
                if target_path.exists():
                    print(f"Skipping {source_path} (duplicate filename).")
                    continue

                # Copy file
                shutil.copy2(source_path, target_path)
                if str(target_path) not in copied_files:
                    copied_files.add(str(target_path))
                    new_files.append(str(target_path))

    # Save updated log
    with open(log_path, "w") as f:
        json.dump(sorted(copied_files), f, indent=2)

    print(f"Copied {len(new_files)} new files. Log updated at {log_file}")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Copy .txt files to target folders")
    parser.add_argument("src", help="Source folder to crawl")
    parser.add_argument("dest", help="Destination folder for normal .txt files")
    parser.add_argument("dest_2d", help="Destination folder for files containing '2D'")
    parser.add_argument("--log", default="copied_files.json", help="Log file to track copied files")
    args = parser.parse_args()

    copy_txt_files(args.src, args.dest, args.dest_2d, args.log)
