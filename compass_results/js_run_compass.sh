#!/bin/bash

#### 80GB mem
#SBATCH --mem=81920 --cpus-per-task=16
#SBATCH --time=24:00:00

source activate compass

dir=/Genomics/pritykinlab/zzhao/metabolic_analysis/metabolic_analysis/data/th_data

compass --data $dir/GSE162300_DFMO_RNA_TPMs.tsv --num-processes 12 --species mus_musculus
