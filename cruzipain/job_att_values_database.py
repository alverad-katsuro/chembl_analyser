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
con = sqlite3.connect("../data/propriedades.db")

# %%
molregno_chembl_id = pd.read_sql_query("select * from dados", con)

# %%
# %%
ids_com_nan = molregno_chembl_id[['chembl_id', 'canonical_smiles']]

# %%

procs = 20
chama_atualiza_in_sql(ids_com_nan, molregno_chembl_id, "dados", '../data/dados_att', procs, procs)
