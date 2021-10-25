# To add a new cell, type '# %%'
# To add a new markdown cell, type '# %% [markdown]'
# %%
import pandas as pd
import sqlite3
import seaborn as sns
from matplotlib import pyplot as plt


# %%
con = sqlite3.connect("data/chembl_29.db")


# %%
molregno = pd.read_sql_query("select molregno, chembl_id from molecule_dictionary", con)


# %%
#molregno = molregno[["molregno", "chembl_id"]]


# %%
molregno


# %%
cruzain_id = pd.read_csv('data/cruzain_id_chemb.csv')


# %%
cruzain_id


# %%
cruzain_id = cruzain_id["Molecule ChEMBL ID"]
cruzain_id


# %%
molregno_chembl_id = molregno.query("chembl_id in @cruzain_id")


# %%
molregno = molregno_chembl_id["molregno"]
molregno


# %%
molregno_chembl_id.head()


# %%
chemb_data = pd.read_sql_query("select * from compound_properties", con)


# %%
chemb_data.duplicated().unique()


# %%
chemb_data.head()


# %%
chemb_data = chemb_data.query("molregno in @molregno_chembl_id.molregno")


# %%
chemb_data_cruzain_inib = pd.merge(molregno_chembl_id, chemb_data, on='molregno')


# %%
pd.set_option('display.max_columns', None)
chemb_data_cruzain_inib


# %%
chemb_data_cruzain_inib.keys()


# %%
chemb_data_cruzain_inib.num_lipinski_ro5_violations.unique()


# %%
chemb_data_cruzain_inib.num_ro5_violations.unique()


# %%
chemb_data_cruzain_inib[chemb_data_cruzain_inib.num_lipinski_ro5_violations == chemb_data_cruzain_inib.num_ro5_violations]


# %%
#chemb_data_cruzain_inib.to_csv("teste.csv")
chemb_data_cruzain_inib['ro3_pass'] = chemb_data_cruzain_inib['ro3_pass'].replace({"Y": 'True', "N": 'False'})


# %%
chemb_data_cruzain_inib.query("num_lipinski_ro5_violations == 0")


# %%
chemb_data_cruzain_inib.query("num_lipinski_ro5_violations == 0 & ro3_pass == 'True'")


# %%
pd.get_dummies(chemb_data_cruzain_inib.query("num_lipinski_ro5_violations == 0 & ro3_pass == 'True'"), columns=["molecular_species"])


# %%


