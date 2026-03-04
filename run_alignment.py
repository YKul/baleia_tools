#########################
# MANIPULATING ALIGNMENTS
#########################

import os
from alignment import extract_seqs_id
from Bio import SeqIO
from Bio import AlignIO

# Run folders for alignment extraction
RUN_FOLDER = "/home/chordata/Desktop/Bielawski_Lab/MSc/results/maleXfemale/genes-multispecies"
RESULTS_FOLDER = "/home/chordata/Desktop/Bielawski_Lab/MSc/results/maleXfemale/alignments"
os.makedirs(RESULTS_FOLDER,exist_ok=True)
LIST_PATH = "/home/chordata/Desktop/Bielawski_Lab/MSc/results/maleXfemale/alignments/maleXfemale.list"


######################################
# Extracting sequences from fasta file
######################################

for sequence_file in os.listdir(RUN_FOLDER):
    if sequence_file.endswith(".fasta"):
        file_path = os.path.join(RUN_FOLDER, sequence_file)
        """
        Yuri's Note:
        Your arguments to extract_seqs_id are file_path and LIST_PATH
        The header for extract_seqs_id expects parameters
        for list_path, seq_folder.

        You have LIST_PATH as a string pointing to a specific .list file,
        but in the method extract_seqs_id you iterate through it
        looking for files ending in '.fasta' as if it were a directory.

        If "small-dataset2.list" is a directory, then '.' can be confusing
        for directory names because it can have special meanings in some contexts.
        Where it appears at the end of a path, and the variable name is LIST_PATH
        instead of LIST_FOLER, it is indistinguishable from a filename'
        """
        extracted_seqs = extract_seqs_id.extract_seqs_id(
            file_path, LIST_PATH)
        with open(f"{RESULTS_FOLDER}/{sequence_file}", "w", encoding="utf-8") as output_file:
            SeqIO.write(extracted_seqs, output_file, "fasta")


######################################
# Extracting sequences from alignment
######################################

# for alignment in os.listdir(RUN_FOLDER):
    # if alignment.endswith(".fasta"):
        # alignment_path = os.path.join(RUN_FOLDER, alignment)

        # extracted_seqs = ExtractSeqsID.extract_seqs_id(
            # alignment_path, LIST_PATH)
        # with open(f"{RESULTS_FOLDER}/{alignment}", "w", encoding="utf-8") as output_file:
            # AlignIO.write(extracted_seqs, output_file, "fasta")
