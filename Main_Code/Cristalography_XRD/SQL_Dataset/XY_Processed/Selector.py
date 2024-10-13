import os
import os
import shutil
from collections import Counter

def copy_frequent_files(min_occurrences, target_dir_name):
    # Get the current directory as the source
    source_dir = os.getcwd()
    
    # Create the target directory relative to the script's location
    script_dir = os.path.dirname(os.path.abspath(__file__))
    target_dir = os.path.join(script_dir, target_dir_name)
    
    # Ensure the target directory exists
    os.makedirs(target_dir, exist_ok=True)

    # Get all file names in the source directory
    all_files = [f for f in os.listdir(source_dir) if f.endswith('.txt')]

    # Count occurrences of each file prefix (everything before the first underscore)
    file_prefixes = [f.split('_')[0] for f in all_files]
    prefix_counts = Counter(file_prefixes)

    # Copy files with prefix occurrences >= min_occurrences
    for file in all_files:
        prefix = file.split('_')[0]
        if prefix_counts[prefix] >= min_occurrences:
            source_path = os.path.join(source_dir, file)
            target_path = os.path.join(target_dir, file)
            shutil.copy2(source_path, target_path)
            print(f"Copied {file} to {target_dir}")

    print(f"Finished copying files with {min_occurrences} or more occurrences.")

# Usage
minimum_occurrences = 6  # Change this to your desired minimum number of occurrences
target_directory_name = "Selected"  # Name of the directory to copy files to

copy_frequent_files(minimum_occurrences, target_directory_name)
