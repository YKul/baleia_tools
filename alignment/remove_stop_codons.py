from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
import sys
from os import listdir
import os

def remove_stop_codons_alignio(input_file, output_file, type="fasta"):
    stop_codons = {"TAA", "TAG", "TGA", "taa", "tag", "tga"}
    
    alignment = AlignIO.read(input_file, type)
    new_records = []
    error_score = 0
    flag_terminal_stop = 0
    for record in alignment:
        seq = str(record.seq)
        
        # Check if sequence length is a multiple of 3
        if len(seq) % 3 != 0:
            print(f"Error: sequence length ({len(seq)}) in ({record.id}) is not a multiple of 3! May contain incomplete codons.")
        else: seq = seq[:-1]
            #return
        
        # Build new sequence without stop codons
        cleaned_seq = []
        for codon_index in range(0, len(seq), 3):
            codon = seq[codon_index:codon_index+3]
            #If a STOP codon is removed from one seq record but not another then 
            #MultipleSeqAlignment throws an errors because of different sequence lengths
            # ^ This will no longer happen as the STOP codon is replaced by "***"
            
            # error_score = 1 STOP in some sequence terminals but not others
            # error_score = 2 STOP detected in reading frame, not terminal codon
            # error_score = 3 More than one STOP detected in a sequence
            if codon in stop_codons:
                #Replaces STOP codons with ***
                cleaned_seq.append("***")
                #get error score for alignment to determine which should be manually checked
                if codon_index+3 == len(seq)-1:
                    #Flags stop codon at end of sequence
                    flag_terminal_stop = 1
                else:
                    if error_score == 2:
                        error_score = 3
                    else:
                        error_score = 2
            else:
                #If some sequences end with STOP codon but not all
                if codon_index+3 == len(seq)-1 and flag_terminal_stop == 1 and error_score == 0:
                    error_score = 1
                cleaned_seq.append(codon)  # Add valid codons to the list
        
        # Convert list to string and create a new sequence record
        new_seq = Seq("".join(cleaned_seq))
        new_records.append(SeqRecord(new_seq, id=record.id, description=""))

    # Write modified alignment
    new_alignment = MultipleSeqAlignment(new_records)
    AlignIO.write(new_alignment, output_file, type)
    
    with open("seq_check_table.csv", 'a') as check_file:
        data = input_file.split("/")[-1] + "," + str(error_score) + "\n"
        check_file.write(data)
    
for sex in ["male","female"]:
    INPUT_DIR = "/home/chordata/Desktop/Bielawski_Lab/MSc/results/"+sex+"/gaptrimmed_alignments/"
    OUTPUT_DIR = "/home/chordata/Desktop/Bielawski_Lab/MSc/results/"+sex+"/finished_alignments/"
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    for file in listdir(INPUT_DIR):
        if not "time" in file:
            # Run the function
            INPUT_FILE = INPUT_DIR+file
            OUTPUT_FILE = OUTPUT_DIR+file
            remove_stop_codons_alignio(INPUT_FILE, OUTPUT_FILE)
