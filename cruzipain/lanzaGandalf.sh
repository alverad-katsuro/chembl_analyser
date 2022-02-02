#!/bin/bash
#SBATCH -J Atualiza_data_base
#SBATCH --time 5-00:00:00
#SBATCH -p cpu               #Fila (partition) a ser utilizada
#SBATCH --ntasks 20

date; pwd;

module load python/ml_39-alfredo

python3 job_att_values_database.py >> arquivo_log.txt

date