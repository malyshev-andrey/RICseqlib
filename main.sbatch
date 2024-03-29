#!/bin/bash
#SBATCH --partition=all
#SBATCH --job-name=RICseqlib
#SBATCH --error=./errors.txt
#SBATCH --output=./logging.txt
#SBATCH --time=3-23:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=300G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<your email>

ENVIRONMENT="./env/snakemake"
SNAKEMAKE_VERSION=8.5.2
rm -rf $ENVIRONMENT

source ~/.bashrc
conda create --yes --prefix $ENVIRONMENT snakemake==$SNAKEMAKE_VERSION
conda activate $ENVIRONMENT

echo 'Snakemake path:'
which snakemake
echo 'Snakemake version:'
snakemake --version

echo "Start: $(date)"
snakemake --profile profile --cores 12
echo "Exit code: ${?}"
echo "Finish: $(date)"
