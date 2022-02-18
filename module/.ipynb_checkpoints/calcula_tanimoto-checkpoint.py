import pandas as pd
import numpy as np
from rdkit.Chem import RDKFingerprint, MolFromSmiles
from rdkit.DataStructs import FingerprintSimilarity
import threading
import sqlite3
from uuid import uuid4
from os import remove

def calcula_tanimoto_module(mol, dataframe, uid):
  columns_name = dataframe.columns.to_list()
  index_smiles = 0
  for colum_rot in columns_name:
    if (colum_rot == 'canonical_smiles'):
        break
    else:
        index_smiles += 1
  dataframe = dataframe.to_numpy().tolist()
  data_fim = []
  while (len(dataframe) != 0):
    dados = dataframe.pop()
    try:
      dados.append(float(f"{FingerprintSimilarity(mol, RDKFingerprint(MolFromSmiles(dados[0]))):.2f}"))   
    except:
      print(f"Error: {dados[0]}\n")
      dados.append(0)
    finally:
      data_fim.append(dados)
  con = sqlite3.connect(uid)
  columns_name.append('i_tanimoto')
  dataframe = pd.DataFrame(data_fim, columns=columns_name)
  dataframe.to_sql("dados", con, if_exists='append', index=False)
  con.close()



def calcula_tanimoto(smiles, dataframe, threads_num) -> pd:
  alvo = RDKFingerprint(MolFromSmiles(smiles))
  jobs = np.array_split(dataframe, threads_num)
  processos = []
  name_data = str(uuid4())[:8]
  for i in range(1, threads_num):
    processThread = threading.Thread(target=calcula_tanimoto_module, args=(alvo, jobs.pop().copy(), f"{name_data}_{i}.db"))
    processThread.start()
    processos.append(processThread)
  calcula_tanimoto_module(alvo, jobs.pop(), f"{name_data}_{0}.db")
  for proc in processos:
    proc.join()
  con = sqlite3.connect(f"{name_data}_{0}.db")
  data_fim = pd.read_sql("select * from dados", con)
  remove(f"{name_data}_{0}.db")
  for i in range(1, threads_num):
    con = sqlite3.connect(f"{name_data}_{i}.db")
    data_fim = pd.concat([data_fim, pd.read_sql("select * from dados", con)])
    remove(f"{name_data}_{i}.db")
  return data_fim