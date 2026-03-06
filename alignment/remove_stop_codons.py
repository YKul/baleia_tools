from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
import sys
from os import listdir
import os
import pandas as pd

def remove_stop_codons_alignio(input_file, output_file, type="fasta"):
    stop_codons = {"TAA", "TAG", "TGA", "taa", "tag", "tga"}
    
    alignment = AlignIO.read(input_file, type)
    new_records = []
    error_score = 0
    flag_terminal_stop = 0

    alignment_summary = {
        'filename':[],
        'record_id':[],
        'check_score':[],
        'num_stops':[],
        'terminal_stop':[],
        'seq_div_3':[]
    }
    
    for record in alignment:
        alignment_summary['filename'].append(input_file.split("/")[-1])
        alignment_summary['record_id'].append(record.id)
    
        seq = str(record.seq)
        
        # Check if sequence length is a multiple of 3
        if len(seq) % 3 != 0:
            print(f"Error: sequence length ({len(seq)}) in ({record.id}) is not a multiple of 3! May contain incomplete codons.")
            alignment_summary['seq_div_3'].append(False)
            alignment_summary['check_score'].append(999)
        else:
            alignment_summary['seq_div_3'].append(True)
            
            # Build new sequence without stop codons
            new_seq = []

            # Build new sequence with stop codons marked
            seq_with_stops = []

            codons = [seq[i:i+3] for i in range(0, len(seq), 3)]
            #codons = set(codons)
            terminal_found = False
            terminal_is_stop = False
            stop_codons_found = 0
            for codon in reversed(codons):
                if not codon == '---':
                    if not terminal_found:
                        terminal_found = True
                        if codon in stop_codons:
                            terminal_is_stop = True
                    if codon in stop_codons:
                        stop_codons_found += 1
                        codon = '***'
                if not codon == '***':        
                    new_seq.append(codon)
                seq_with_stops.append(codon)
            alignment_summary['num_stops'].append(stop_codons_found)
            alignment_summary['terminal_stop'].append(terminal_is_stop)
            
            # check_score = 0 One STOP in some terminal only
            # check_score = 1 No STOP codons
            # check_score = 2 STOP detected in reading frame, not terminal codon
            # check_score = 3 More than one STOP detected in a sequence

            if stop_codons_found > 1:
                alignment_summary['check_score'].append(3)
            elif stop_codons_found == 0:
                alignment_summary['check_score'].append(1)
            else:
                if terminal_is_stop:
                    alignment_summary['check_score'].append(0)
                else:
                    alignment_summary['check_score'].append(2)
            
            
        # Convert list to string and create a new sequence record
        new_seq = Seq("".join(reversed(new_seq)))
        #new_seq = Seq("".join(reversed(seq_with_stops)))
        new_records.append(SeqRecord(new_seq, id=record.id, description=""))

    # Write modified alignment
    new_alignment = MultipleSeqAlignment(new_records)
    AlignIO.write(new_alignment, output_file, type)
    
    pd.DataFrame.from_dict(alignment_summary).to_csv(alignment_summary['filename'][0].split('.')[0]+".csv")


#This is implemented only for my own local testing - Yuri
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
                INPUT_DIR = line.split('=')[-1]+sex+"/gaptrimmed_alignments/"
                OUTPUT_DIR = line.split('=')[-1]+sex+"/finished_alignments/"
                #print(INPUT_DIR)
                #print(OUTPUT_DIR)
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    for file in listdir(INPUT_DIR):
        #skips the files storing runtimes
        if not "time" in file and not os.path.isdir(INPUT_DIR+file):
            # Run the function
            INPUT_FILE = INPUT_DIR+file
            OUTPUT_FILE = OUTPUT_DIR+file
            remove_stop_codons_alignio(INPUT_FILE, OUTPUT_FILE)
