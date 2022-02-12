#!/bin/bash
#SBATCH -J Busca_smiles_bug
#SBATCH --time 5-00:00:00
#SBATCH -p cpu               #Fila (partition) a ser utilizada
#SBATCH --ntasks 20

date; pwd;

module load python/ml_39-alfredo

#python3 job_att_values_database.py
python3 job_02.py

date
