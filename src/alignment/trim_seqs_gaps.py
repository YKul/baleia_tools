from Bio import AlignIO
from Bio import SeqIO
import os
from os import listdir

for sex in ["male","female"]:
    with open ("../config.txt",'r') as config_file:
        for line in config_file:
            #Sanitize
            line = line.strip().replace(" ","")
            if "results_dir" in line:                
                #Adds the trailing / in case it's forgotten in the config.txt
                #to just avoid errors altogether
                if line[-1] != '/':
                    line += '/'
                INPUT_DIR = line.split('=')[-1]+sex+"/taper_corrected_alignments/"
                OUTPUT_DIR = line.split('=')[-1]+sex+"/gaptrimmed_alignments/"
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    MAX_GAP_PERCENTAGE = 50.0  # maximum allowed gap percentage
    for file in listdir(INPUT_DIR):
        print(INPUT_DIR+file)
        if not "time" in file:
            # Read the alignment
            alignment = AlignIO.read(INPUT_DIR+file, "fasta")
            alignment_length = alignment.get_alignment_length()
            max_gaps = (MAX_GAP_PERCENTAGE / 100.0) * alignment_length

            # Filter sequences
            filtered_sequences = []
            for record in alignment:
                if record.seq.count("-") <= max_gaps:
                    filtered_sequences.append(record)

            # Write the filtered alignment
            SeqIO.write(filtered_sequences, OUTPUT_DIR+file, "fasta")


