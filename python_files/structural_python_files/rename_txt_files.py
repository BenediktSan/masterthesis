import os

def rename_txt_files():
    # Get the current working directory
    current_dir = os.getcwd()
    
    # Iterate through each file in the directory
    for filename in os.listdir(current_dir):
        if filename.endswith(".txt") and not filename.endswith("_old.txt"):
            base_name = filename[:-4]  # Remove .txt
            new_name = f"{base_name}_old.txt"
            os.rename(filename, new_name)
            print(f"Renamed: {filename} → {new_name}")

if __name__ == "__main__":
    rename_txt_files()
