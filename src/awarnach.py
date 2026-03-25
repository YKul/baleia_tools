import os
from src import awarnach_config
from posixpath import join as posixjoin

cfg = awarnach_config.Config()
RESULTS_DIRS= cfg.getResultsDirs()
MACSE_PATH = cfg.getMACSE()
JAVA_PATH = cfg.getJAVA()

def _createHeader(array_size=0,job_name="unnamed_job",stream_output="",error_output="",i=0):
    header = f"\
#! /bin/bash\n\
#Path to shell\n\
#$ -S /bin/bash\n\
. /etc/profile\n\
#Output directory\n\
#$ -o {cfg.getStreamOutput(i)}\n\
#$ -e {cfg.getErrorOutput(i)}\n\
#e-mail\n\
#$ -M {cfg.getEmail()}\n\
#Job Name\n\
#$ -N {job_name}\n\
#Specify queue\n\
#$ -q 256G-batch@awarnach*\n"
    if array_size > 0:
        header += f"#Array Size\n#$ -t 1-{array_size}\n"
    header += "\n"

    return header

#TODO: verify this still works
def makeJobMAFFT(multispecies_gene_directory,cluster_jobs):
    for i in range(len(RESULTS_DIRS)):
        array_size = len(os.listdir(multispecies_gene_directory[i]))
        header = _createHeader(array_size,"align_MAFFT"+str(i),i)
        body = f"\
dataset=$((SGE_TASK_ID-1))\n\
# Uses the gene names in the genes-multispecies directory to avoid added suffixes\
input_dir='{posixjoin(RESULTS_DIRS[i],"genes-multispecies")}'\n\
input_files=({posixjoin(RESULTS_DIRS[i],"genes-multispecies",'*.fasta')})\n\
file_handle=$(basename ${{input_files[$dataset]}})\n\
no_extension=${{file_handle%.*}}\n\
\n\
input_dir={posixjoin(RESULTS_DIRS[i],"macse_aligned",'$no_extension')}\n\
input_file={posixjoin(RESULTS_DIRS[i],"macse_aligned",'$no_extension','$no_extension'+"'_nuc_AA.fasta'")}\n\
output_dir={posixjoin(RESULTS_DIRS[i],"MAFFT_aligned",'$no_extension')}\n\
\n\
if [ ! -d $output_dir ]; then\n\
    mkdir -p $output_dir\n\
fi\n\n"

        body += "/usr/bin/time -v mafft --globalpair --maxiterate 1000 $input_dir/$file_handle > $output_dir/$file_handle 2> $output_dir/$file_handle.time.global"

        with open(os.path.join(cluster_jobs[i],"nuc_MAFFT-group"+str(i)+".sh"),'w') as output:
            output.write(header+body)

def makeJobTrimNonHomoFrags(multispecies_gene_directory,cluster_jobs):
    for i in range(len(RESULTS_DIRS)):
        array_size = len(os.listdir(multispecies_gene_directory[i]))
        header = _createHeader(array_size,"TNH"+str(i),i)
        body = f"\
dataset=$((SGE_TASK_ID-1))\n\
input_dir='{posixjoin(RESULTS_DIRS[i],"genes-multispecies")}'\n\
input_files=({posixjoin(RESULTS_DIRS[i],"genes-multispecies",'*.fasta')})\n\
file_handle=$(basename ${{input_files[$dataset]}})\n\
no_extension=${{file_handle%.*}}\n\
output_dir={posixjoin(RESULTS_DIRS[i],"trimmed_nonhomo",'$no_extension')}\n\
\n\
if [ ! -d $output_dir ]; then\n\
    mkdir -p $output_dir\n\
fi\n\n\
\n\n"

        body += f"/usr/bin/time -v {JAVA_PATH} -jar '{MACSE_PATH}' -prog trimNonHomologousFragments\
 -seq $input_dir/$file_handle -out_AA $output_dir/$no_extension'_AA.fasta' -out_NT $output_dir/$no_extension'_NT.fasta'\
 -out_trim_info $output_dir/$no_extension'_homol_fiter.csv' -out_mask_detail\
 $output_dir/$no_extension'_NonHomolFilter_NT_mask_detail.fasta' 2> $output_dir/$file_handle.time.global"
        with open(os.path.join(cluster_jobs[i],"nuc_trimnonhomo-group"+str(i)+".sh"),'w') as output:
            output.write(header+body)    

#Remember you're specifying directories here because you need to make a POSIX path
def makeJobMACSE(multispecies_gene_directory,cluster_jobs):
    for i in range(len(RESULTS_DIRS)):
        array_size = len(os.listdir(multispecies_gene_directory[i]))
        header = _createHeader(array_size,"macse"+str(i),i)
        body = f"\
dataset=$((SGE_TASK_ID-1))\n\
# Uses the gene names in the genes-multispecies directory to avoid added suffixes\
input_dir='{posixjoin(RESULTS_DIRS[i],"genes-multispecies")}'\n\
input_files=({posixjoin(RESULTS_DIRS[i],"genes-multispecies",'*.fasta')})\n\
file_handle=$(basename ${{input_files[$dataset]}})\n\
no_extension=${{file_handle%.*}}\n\
\n\
input_dir={posixjoin(RESULTS_DIRS[i],"trimmed_nonhomo",'$no_extension')}\n\
input_file={posixjoin(RESULTS_DIRS[i],"trimmed_nonhomo",'$no_extension','$no_extension'+"'_NT.fasta'")}\n\
output_dir={posixjoin(RESULTS_DIRS[i],"macse_aligned",'$no_extension')}\n\
\n\
if [ ! -d $output_dir ]; then\n\
    mkdir -p $output_dir\n\
fi\n\n"

        body += f"/usr/bin/time -v {JAVA_PATH} -jar '{MACSE_PATH}' -prog alignSequences\
 -seq $input_file -out_NT $output_dir/$no_extension'_NT.fasta' -out_AA $output_dir/$no_extension'_AA.fasta'\
  2> $output_dir/$file_handle.time.global\n"


        with open(os.path.join(cluster_jobs[i],"nuc_macse_align"+str(i)+".sh"),'w') as output:
            output.write(header+body)    