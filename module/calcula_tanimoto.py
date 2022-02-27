import pandas as pd
import numpy as np
from rdkit.Chem import RDKFingerprint, MolFromSmiles, MolToSmiles, SanitizeMol
from rdkit.DataStructs import FingerprintSimilarity
from multiprocessing import Process
import sqlite3
from uuid import uuid4
from os import remove, mkdir, rmdir

import sqlite3
from uuid import uuid4
from os import remove

def sanitize(mol):
  try:
      SanitizeMol(mol)
      return mol
  except Exception as e:
      print(f"{MolToSmiles(mol)} {e}")

def calcula_tanimoto_module(mol, columns_name, numpy_arr, uid):
  if ("canonical_smiles" in columns_name):
    smiles_index = columns_name.index("canonical_smiles")
  else:
    smiles_index = columns_name.index("Smiles")
  data_fim = []
  numpy_arr = numpy_arr.tolist()
  while (len(numpy_arr) != 0):
    dados = numpy_arr.pop()
    try:
      dados.append(float(f"{FingerprintSimilarity(mol, RDKFingerprint(sanitize(MolFromSmiles(dados[smiles_index])))):.2f}"))  
    except:
      print(f"Error: {dados[smiles_index]}\n")
      dados.append(0)
    finally:
      data_fim.append(dados)
  con = sqlite3.connect(uid)
  columns_name.append('i_tanimoto')
  dataframe = pd.DataFrame(data_fim, columns=columns_name)
  dataframe.to_sql("dados", con, if_exists='append', index=False)
  con.close()



def calcula_tanimoto(smiles, dataframe, threads_num) -> pd:
  mol = RDKFingerprint(MolFromSmiles(smiles))
  columns_name = list(dataframe.columns)
  dataframe = np.array_split(dataframe.to_numpy().tolist(), threads_num)
  processos = []
  name_data = str(uuid4())[:8]
  mkdir(name_data)
  for i in range(threads_num):
    p = Process(target=calcula_tanimoto_module, args=(mol, columns_name, dataframe.pop(), f"{name_data}/{name_data}_{i}.db"))
    p.start()
    processos.append(p)
  for proc in processos:
    proc.join()
  con = sqlite3.connect(f"{name_data}/{name_data}_{0}.db")
  data_fim = pd.read_sql("select * from dados", con)
  remove(f"{name_data}/{name_data}_{0}.db")
  for i in range(1, threads_num):
    con = sqlite3.connect(f"{name_data}/{name_data}_{i}.db")
    data_fim = pd.concat([data_fim, pd.read_sql("select * from dados", con)])
    remove(f"{name_data}/{name_data}_{i}.db")
  rmdir(name_data)
  return data_fim