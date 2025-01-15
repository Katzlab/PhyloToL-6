'''
Author & Date: Adri K. Grow + ChatGPT, Nov 11th 2024
Motivation: assess and rename assembled transcript files for PTL6p1 (or EPv1p1)
Intention: warn if any 'transcripts.fasta' files are missing or empty for an LKH, otherwise rename and copy them with their newly assigned 10-digit code by LKH
Input:
    - a base directory containing subdirectories for each LKH assembled transcriptome, named 'WTA_LKH<xxxx>', each containing a 'transcripts.fasta' file
    - a mapping .txt file with LKH#s tab-separated with corresponding 10-digit codes
Output: 
    - a folder named 'renamed_transcripts' with 'transcripts.fasta' files now named by 10-digit codes; e.g. "Sr_rh_Ro04_assembledTranscripts.fasta"
Dependencies: python3
Usage: python3 ProcessAndRenamePlateTranscripts_AKG.py <assembled transcriptomes directory> <mapping_file.txt>
'''

import os
import shutil
import sys

def read_lkh_mapping(mapping_file):
    """Reads the LKH number to 10-digit code mapping from a file."""
    mapping = {}
    with open(mapping_file, 'r') as file:
        for line in file:
            lkh_number, code = line.strip().split('\t')
            mapping[lkh_number] = code
    return mapping

def process_directory(base_dir, mapping, output_dir):
    """Iterates over all subdirectories in base_dir, processes transcripts.fasta files."""
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)  # Create output directory if it doesn't exist
    
    # Iterate over subdirectories in base_dir
    for folder_name in os.listdir(base_dir):
        folder_path = os.path.join(base_dir, folder_name)

        if os.path.isdir(folder_path) and folder_name.startswith("WTA_LKH"):
            lkh_number = folder_name.split('_')[1]  # Extract LKH number from folder name
            transcripts_file = os.path.join(folder_path, 'transcripts.fasta')
            
            # Check if the file exists and is not empty
            if not os.path.isfile(transcripts_file):
                print(f"    WARNING: file 'transcripts.fasta' is missing in folder {folder_name}.")
                continue  # Skip to the next folder
            
            if os.path.getsize(transcripts_file) == 0:
                print(f"    WARNING: file 'transcripts.fasta' is empty in folder {folder_name}.")
                continue  # Skip to the next folder
            
            # If file exists and is not empty, proceed with renaming and copying
            if lkh_number in mapping:
                new_name = f"{mapping[lkh_number]}_assembledTranscripts.fasta"
                output_path = os.path.join(output_dir, new_name)
                
                # Copy the file with the new name
                shutil.copy(transcripts_file, output_path)
                #print(f"Copied and renamed {transcripts_file} to {output_path}.")
            else:
                print(f"Notification: No 10-digit code found for LKH number {lkh_number} in folder {folder_name}.")

def main():
    # Check if the user provided the correct number of command-line arguments
    if len(sys.argv) != 3:
        print("Usage: python script.py <base_dir> <mapping_file>")
        sys.exit(1)
    
    # Get the base directory and mapping file from command line arguments
    base_dir = sys.argv[1]
    mapping_file = sys.argv[2]
    
    # Check if the directories/files exist
    if not os.path.isdir(base_dir):
        print(f"Error: The directory '{base_dir}' does not exist.")
        sys.exit(1)

    if not os.path.isfile(mapping_file):
        print(f"Error: The file '{mapping_file}' does not exist.")
        sys.exit(1)
    
    # Output directory will be created in the current working directory
    output_dir = os.path.join(os.getcwd(), "renamed_transcripts")
    
    # Read the LKH to code mapping
    mapping = read_lkh_mapping(mapping_file)
    
    # Process the directories
    process_directory(base_dir, mapping, output_dir)

if __name__ == "__main__":
    main()
