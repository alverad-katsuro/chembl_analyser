#!/usr/bin/env python
'''
Chemical data about a molecule.

Molecules are defined by SMILES strings. Can work out logP values, Lipinski's 
rules, etc...

Uses rdkit
'''

from rdkit import Chem
from rdkit.Chem import Crippen
from rdkit.Chem import Lipinski
from rdkit.Chem import Descriptors
from rdkit.Chem.QED import properties


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
    
def lipinski_pass(num_hdonors, num_hacceptors, mol_weight, mol_logp):
    '''
    Wraps around lipinski trial, but returns a simple pass/fail True/False
    '''
    passed, failed = lipinski_trial(num_hdonors, num_hacceptors, mol_weight, mol_logp)
    if failed:
        return {'Pass?': False, "N_RO5": float(len(failed))}
    else:
        return {'Pass?': True, "N_RO5": float(len(failed))}


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
    resultados['aromatic_rings'] = float(f"{Chem.GetSSSR(mol):.2f}")
    resultados['psa'] = float(f"{Chem.MolSurf.TPSA(mol):.2f}")
    

    teste_lin = lipinski_pass(resultados['hbd_lipinski'], resultados['hba_lipinski'], resultados['mw_freebase'], resultados['alogp'])

    resultados['num_lipinski_ro5_violations'] = f"{teste_lin['N_RO5']:.4f}"

    return resultados

def atualiza_data_frame_com_lipinski(dataframe, keys):
  for chemb, smiles in zip(keys.chembl_id, keys.canonical_smiles):
    propriedades = verifica_lipinski(smiles)
    for key in propriedades.keys():
        dataframe.loc[dataframe.chembl_id == chemb, key] = propriedades[key]