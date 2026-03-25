import subprocess
from Bio import SeqIO
from io import StringIO
import os

#Algorithm https://mhmmdbduh.medium.com/sequence-alignment-using-mafft-in-python-335a4a84a6d4
def localMAFFT(multispecies_gene_directory, alignment_output_dir):
    for file in os.listdir(multispecies_gene_directory):
        seq_str = ""
        filepath = os.path.join(multispecies_gene_directory, file)
        fastas = list(SeqIO.parse(filepath,format='fasta'))
        for seq in fastas:
            seq_str += ">"+seq.description+"\n"
            seq_str += str(seq.seq)+"\n"
        print(f"Aligning: {filepath}")
        child = subprocess.Popen(['mafft', '-'], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        child.stdin.write(seq_str.encode())
        child_out = child.communicate()[0].decode('utf8')
        seq_ali = list(SeqIO.parse(StringIO(child_out), 'fasta'))
        child.stdin.close()

        output_dir = os.path.join(alignment_output_dir,file)
        with open(output_dir,'w') as output_handle:
            SeqIO.write(seq_ali,output_handle,'fasta')

