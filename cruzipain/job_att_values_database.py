# %%
# To add a new cell, type '# %%'
# To add a new markdown cell, type '# %% [markdown]'
# %%
import pandas as pd
import sqlite3
import sys, os
sys.path.append(os.path.abspath("../module"))
from lipinski import * 

# %%
pd.options.display.max_columns = 30

# %%
con = sqlite3.connect("../data/chembl_29.db")

# %%
molregno_chembl_id = pd.read_sql_query("SELECT compound_properties.*, molecule_dictionary.chembl_id, compound_structures.canonical_smiles FROM compound_properties INNER JOIN compound_structures ON compound_structures.molregno=compound_properties.molregno INNER JOIN molecule_dictionary ON molecule_dictionary.molregno=compound_properties.molregno;", con)

# %%
molregno_chembl_id = molregno_chembl_id.drop(['ro3_pass', 'hba', 'hbd', 'num_ro5_violations', 'full_molformula'], axis=1)

# %%
molregno_chembl_id = pd.merge(molregno_chembl_id.chembl_id, molregno_chembl_id)

# %%
molregno_chembl_id = molregno_chembl_id.drop(['molregno'], axis=1)

# %%
index_duplicates = molregno_chembl_id[molregno_chembl_id.canonical_smiles.duplicated()].sort_values(by=['canonical_smiles']).index

# %%
molregno_chembl_id = molregno_chembl_id.drop(labels=index_duplicates, axis=0)
molregno_chembl_id

# %%
ids_com_nan = molregno_chembl_id.loc[molregno_chembl_id.isnull().any(axis=1)][['chembl_id', 'canonical_smiles']]

# %%
ids_com_nan

# %%

chama_atualiza_in_sql(ids_com_nan, molregno_chembl_id, "dados", '../data/dados_att', 6*3, 6)
