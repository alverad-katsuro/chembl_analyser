#!/usr/bin/env python
'''
Chemical data about a molecule.

Molecules are defined by SMILES strings. Can work out logP values, Lipinski's 
rules, etc...

Uses rdkit
'''

from xmlrpc.client import Boolean
from rdkit import Chem
from rdkit.Chem import Crippen
from rdkit.Chem import Lipinski
from rdkit.Chem import Descriptors
import numpy as np
from multiprocessing import Process

import queue
import sqlite3

class SmilesError(Exception): pass

def log_partition_coefficient(smiles):
  '''
  Returns the octanol-water partition coefficient given a molecule SMILES 
  string
  '''
  try:
      mol = Chem.MolFromSmiles(smiles)
  except Exception as e:
      print(f"Ocorreu um erro {e}")
      raise SmilesError('%s returns a None molecule' % smiles)
      
  return Crippen.MolLogP(mol)
    
def lipinski_trial(num_hdonors, num_hacceptors, mol_weight, mol_logp):
  '''
  Returns which of Lipinski's rules a molecule has failed, or an empty list
  
  Lipinski's rules are:
  Hydrogen bond donors <= 5
  Hydrogen bond acceptors <= 10
  Molecular weight < 500 daltons
  logP < 5
  '''
  passed = []
  failed = []
  
  if num_hdonors > 5:
      failed.append('Over 5 H-bond donors, found %s' % num_hdonors)
  else:
      passed.append('Found %s H-bond donors' % num_hdonors)
      
  if num_hacceptors > 10:
      failed.append('Over 10 H-bond acceptors, found %s' \
      % num_hacceptors)
  else:
      passed.append('Found %s H-bond acceptors' % num_hacceptors)
      
  if mol_weight >= 500:
      failed.append('Molecular weight over 500, calculated %s'\
      % mol_weight)
  else:
      passed.append('Molecular weight: %s' % mol_weight)
      
  if mol_logp >= 5:
      failed.append('Log partition coefficient over 5, calculated %s' \
      % mol_logp)
  else:
      passed.append('Log partition coefficient: %s' % mol_logp)
  
  return passed, failed
    
def lipinski_pass_modulo(num_hdonors, num_hacceptors, mol_weight, mol_logp):
  '''
  Wraps around lipinski trial, but returns a simple pass/fail True/False
  '''
  passed, failed = lipinski_trial(num_hdonors, num_hacceptors, mol_weight, mol_logp)
  if failed:
      return {'Pass?': False, "N_RO5": float(len(failed))}
  else:
      return {'Pass?': True, "N_RO5": float(len(failed))}

def lipinski_pass(smiles) -> Boolean:
  mol = Chem.MolFromSmiles(smiles)
  if mol is None:
      raise Exception('%s is not a valid SMILES string' % smiles)
  
  resultados = {}
  
  resultados['hbd_lipinski'] = float(f"{Lipinski.NumHDonors(mol):.4f}")
  resultados['hba_lipinski'] = float(f"{Lipinski.NumHAcceptors(mol):.4f}")
  resultados['mw_freebase'] = float(f"{Descriptors.MolWt(mol):.4f}")
  resultados['alogp'] = float(f"{Crippen.MolLogP(mol):.4f}")
  resultados['rtb'] = float(f"{Lipinski.NumRotatableBonds(mol):.2f}")
  resultados['aromatic_rings'] = float(f"{Chem.GetSSSR(mol):.1f}")
  resultados['psa'] = float(f"{Chem.MolSurf.TPSA(mol):.2f}")
  resultados['heavy_atoms'] = float(f"{mol.GetNumHeavyAtoms():.1f}")
  resultados['qed_weighted'] = float(f"{Chem.QED.qed(mol):.2f}")
  return lipinski_pass_modulo(resultados['hbd_lipinski'], resultados['hba_lipinski'], resultados['mw_freebase'], resultados['alogp'])

def lipinski_pass_dataframe(dataframe) -> list:
  columns_name = list(dataframe.columns)
  if ("canonical_smiles" in columns_name):
    col = "canonical_smiles"
  else:
    col = "Smiles"
  dataframe = dataframe[col].to_numpy().tolist()
  valores = []
  for smiles in dataframe:
    try:
      valores.append(lipinski_pass(smiles)["Pass?"])
    except:
      valores.append(np.nan)
  return valores
  

def verifica_lipinski(smiles) -> dict:
  '''
  Returns which of Lipinski's rules a molecule has failed, or an empty list
  
  Lipinski's rules are:
  Hydrogen bond donors <= 5
  Hydrogen bond acceptors <= 10
  Molecular weight < 500 daltons
  logP < 5
  '''
  
  mol = Chem.MolFromSmiles(smiles)
  if mol is None:
      raise Exception('%s is not a valid SMILES string' % smiles)
  
  resultados = {}
  
  resultados['hbd_lipinski'] = float(f"{Lipinski.NumHDonors(mol):.4f}")
  resultados['hba_lipinski'] = float(f"{Lipinski.NumHAcceptors(mol):.4f}")
  resultados['mw_freebase'] = float(f"{Descriptors.MolWt(mol):.4f}")
  resultados['alogp'] = float(f"{Crippen.MolLogP(mol):.4f}")
  resultados['rtb'] = float(f"{Lipinski.NumRotatableBonds(mol):.2f}")
  resultados['aromatic_rings'] = float(f"{Chem.GetSSSR(mol):.1f}")
  resultados['psa'] = float(f"{Chem.MolSurf.TPSA(mol):.2f}")
  resultados['heavy_atoms'] = float(f"{mol.GetNumHeavyAtoms():.1f}")
  resultados['qed_weighted'] = float(f"{Chem.QED.qed(mol):.2f}")

  teste_lin = lipinski_pass(resultados['hbd_lipinski'], resultados['hba_lipinski'], resultados['mw_freebase'], resultados['alogp'])

  resultados['num_lipinski_ro5_violations'] = float(f"{teste_lin['N_RO5']:.2f}")

  return resultados

def modulo_atualiza(jobs, dataframe):
  while True:
    try:
      chemb, smile = jobs.get(timeout=3)  # 3s timeout
      propriedades = verifica_lipinski(smile)
      for key in list(propriedades):
        dataframe.loc[dataframe.chembl_id == chemb, key] = propriedades[key]
    except queue.Empty:
      return
    jobs.task_done()

def counting_threads(jobs, dataframe, threads_num):
  processos = []    
  for _ in range(threads_num):
    processThread = Process(target=modulo_atualiza, args=(jobs, dataframe))
    processThread.start()
    processos.append(processThread)
  for proc in processos:
    proc.join()

def atualiza_data_frame_com_lipinski(ids_com_nan, dataframe, threads_num):
  jobs = queue.Queue()
  for chemb, smiles in zip(ids_com_nan.chembl_id, ids_com_nan.canonical_smiles):
    jobs.put_nowait([chemb, smiles])
  counting_threads(jobs, dataframe.copy(), threads_num)

def modulo_atualiza_in_sql(jobs, dataframe, nome_table, con_dir, i):
  con = sqlite3.connect(f"{con_dir}/dados_atualizados_{i}.db")
  while (len(jobs) != 0):
    chemb, smile = jobs.pop()
    try:
        propriedades = verifica_lipinski(smile)
        dataframe.loc[dataframe.chembl_id == chemb].copy()
        for key in list(propriedades):
          dataframe.loc[dataframe.chembl_id == chemb, key] = propriedades[key]
    except:
        error = open("smile_erros.txt", "a")
        error.write(f"{smile}\n")
        error.close()
  dataframe.to_sql(nome_table, con, if_exists='append', index=False)
  con.close()

def atualiza_data_frame_com_lipinski_in_sql(ids_com_nan, dataframe, nome_table, con, threads_num):
  jobs = np.array_split(ids_com_nan, threads_num)
  processos = []    
  for i in range(1, threads_num):
    job = jobs.pop()
    processThread = Process(target=modulo_atualiza_in_sql, args=(job.to_numpy().tolist(), dataframe.query("chembl_id in @job['chembl_id']").copy(), nome_table, con, i))
    processThread.start()
    processos.append(processThread)
  modulo_atualiza_in_sql(job.to_numpy().tolist(), dataframe.query("chembl_id in @job['chembl_id']").copy(), nome_table, con, 0)
  for proc in processos:
    proc.join()

def chama_atualiza_in_sql(ids_com_nan, dataframe, nome_table, con_dir, quantidades, threads_num):
  dataframe = dataframe.query("chembl_id in @ids_com_nan['chembl_id']")
  log = open("arquivo.log", "a")
  log.write("Começou\n")
  log.close()
  for i in range(quantidades):
    log = open("arquivo.log", "a")
    print(f"Parte {i + 1} de {quantidades}\n")
    log.write(f"Parte {i + 1} de {quantidades}\n")
    log.close()
    tamanho_ini = int(i * len(ids_com_nan)/quantidades)
    tamanho_fim = int((i + 1) * len(ids_com_nan)/quantidades)
    atualiza_data_frame_com_lipinski_in_sql(ids_com_nan.iloc[tamanho_ini:tamanho_fim], dataframe.iloc[tamanho_ini:tamanho_fim], nome_table, con_dir, threads_num)
  log = open("arquivo.log", "a")
  log.write("Terminou\n")
  log.close()
