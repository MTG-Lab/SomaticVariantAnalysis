#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=rumika.mascarenhas@ucalgary.ca
#SBATCH --job-name=job.mutect.filtering
#SBATCH --time=10:00:00
#SBATCH --partition=cpu2021
#SBATCH --mem-per-cpu=9500
#SBATCH --output=job.out
#SBATCH --array=10
#SBATCH --error=job.err
#SBATCH --mail-type=END


# Specify the path to the config file
config=config_samples.txt


# Extract the sampleID for the current $SLURM_ARRAY_TASK_ID
sampleid=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)

# Extract the affected_sample for the current $SLURM_ARRAY_TASK_ID
affected=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $config)

# Extract the unaffected_sample for the current $SLURM_ARRAY_TASK_ID
unaffected=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $4}' $config)

# Extract the filename for the current $SLURM_ARRAY_TASK_ID
filename=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $5}' $config)

python VariantFiltration.py ${sampleid} ${affected} ${unaffected}

cd /work/mtgraovac_lab/tools/bcftools-1.10.2/

./bcftools sort /work/mtgraovac_lab/Rumika/SeqRound2/MutectRuns/${filename} -o /work/mtgraovac_lab/Rumika/SeqRound2/MutectRuns/sortedVCF/${filename}

