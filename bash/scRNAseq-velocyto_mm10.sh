#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=velocyto_%A_%a.txt
#SBATCH --job-name=velocyto
#SBATCH --mem=256G  # memory requested, units available: K,M,G,T

conda activate r4.1.1
spack load -r samtools@1.8
samtools --version
echo $SLURM_ARRAY_TASK_ID
#---------------------Variables to be set-------------------------#
PROJECT_NAME="scRNAseq-GBM"
path=/athena/elementolab/scratch/yah2014/Projects/${PROJECT_NAME}/data/bam
file_folder=$(ls ${path} | tail -n +${SLURM_ARRAY_TASK_ID}| head -1) # Uses job array for each sample in the folder
rmsk_gtf=/athena/elementolab/scratch/yah2014/Indexed_genome/mm10_rmsk.gtf
genes_gtf=/athena/elementolab/scratch/yah2014/Indexed_genome/refdata-gex-mm10-2020-A/genes/genes.gtf
echo "path="
echo "$path"
echo " "
echo $(ls -l $path/$file_folder/possorted_genome_bam.bam)
echo $(ls -l $rmsk_gtf)
echo $(ls -l $genes_gtf)

# don't sort bam file separately form Velocyto
# ==== check if loom exist ======
if test -f "${path%bam}velocyto/$file_folder.loom"
then
    return 0;
    echo "loom file exist"
else
    echo "loom file does not exist. Continue velocyto run10x."
fi

#----------------successfully File-------------------
echo "to sort by cellID"
echo $(ls -l $path/$file_folder/possorted_genome_bam.bam)
echo Pipestance completed successfully! > $path/$file_folder/_log

#-----------velocyto Command--------------------------------#
echo "Processing velocyto run10x"
echo " "
echo "-------------------------------- "
echo "Processing $file_folder"
cd $path/$file_folder
velocyto run -b $path/$file_folder/barcodes.tsv -e $file_folder -m $rmsk_gtf $path/$file_folder/possorted_genome_bam.bam $genes_gtf

# ==== check if loom exist ======
echo "velocyto output files:"
echo " "
echo $(ls -l $path/$file_folder/velocyto/$file_folder.loom)
if test -f "$path/$file_folder/velocyto/$file_folder.loom"
then
    echo "loom file exists"
    rsync -rav $path/$file_folder/velocyto/$file_folder.loom ${path%bam}velocyto

    echo "velocyto run10x Complished"
else
    echo "loom file does not exist."
    echo " "
    echo "velocyto run10x Failed."
fi
