#!/bin/bash
#SBATCH --job-name=
#SBATCH --partition=
#SBATCH --output=%x.out
#SBATCH --error=%x.err
#SBATCH --mem=
#SBATCH --cpus-per-task=


# Set the paths to your input directory and rankedlineage.dmp file
INPUT_DIR=
RANKEDLINEAGE=
OUTPUT=


# Run the Python script
srun python 3_local_bold_processing.py $INPUT_DIR $RANKEDLINEAGE $OUTPUT

echo complete!
