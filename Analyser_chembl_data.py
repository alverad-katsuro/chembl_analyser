# To add a new cell, type '# %%'
# To add a new markdown cell, type '# %% [markdown]'
# %%
import pandas as pd
import sqlite3
import seaborn as sns
from matplotlib import pyplot as plt
from sklearn.cluster import KMeans
from sklearn import preprocessing
import numpy as np
from matplotlib import pyplot as ptl
import seaborn as sns


# %%
con = sqlite3.connect("data/chembl_29.db")


# %%
molregno_chembl_id = pd.read_sql_query("SELECT compound_properties.*, molecule_dictionary.chembl_id, compound_structures.canonical_smiles FROM compound_properties INNER JOIN compound_structures ON compound_structures.molregno=compound_properties.molregno INNER JOIN molecule_dictionary ON molecule_dictionary.molregno=compound_properties.molregno;", con)


# %%
molregno_chembl_id


# %%
cruzain_id = pd.read_csv('data/cruzain_id_chemb.csv')


# %%
cruzain_id


# %%
cruzain_id = cruzain_id["Molecule ChEMBL ID"]
cruzain_id


# %%
molregno_chembl_id = molregno_chembl_id.query("chembl_id in @cruzain_id")


# %%
molregno_chembl_id.head()


# %%
pd.set_option('display.max_columns', None)


# %%
molregno_chembl_id.num_lipinski_ro5_violations.unique()


# %%
molregno_chembl_id.num_ro5_violations.unique()


# %%
molregno_chembl_id[molregno_chembl_id.num_lipinski_ro5_violations == molregno_chembl_id.num_ro5_violations]


# %%
#chemb_data_cruzain_inib.to_csv("teste.csv")
molregno_chembl_id['ro3_pass'] = molregno_chembl_id['ro3_pass'].replace({"Y": 'True', "N": 'False'})


# %%
molregno_chembl_id.query("num_lipinski_ro5_violations == 0")


# %%
molregno_chembl_id.query("num_lipinski_ro5_violations == 0 & ro3_pass == 'True'")


# %%
prep_clust = pd.get_dummies(molregno_chembl_id.query("num_lipinski_ro5_violations == 0 & ro3_pass == 'True'"), columns=["molecular_species"])
tab_usada = prep_clust
prep_clust = prep_clust.drop(['chembl_id','molregno', "num_ro5_violations", "ro3_pass", "full_molformula", "canonical_smiles"], axis=1)
prep_clust.head()

# %%
prep_clust = prep_clust.replace({np.nan: "7"})
prep_clust.head()

#%%
x = preprocessing.normalize(prep_clust, norm='l1')
cluster = KMeans(n_clusters=3, random_state=5)
cluster.fit(x).labels_

#%%
plt.figure(figsize=(10,7))
plt.scatter(x=prep_clust[['mw_freebase']], y=prep_clust.psa, c=cluster.fit(x).labels_, s=50)
plt.xlabel('mw_freebase', fontsize=18)
plt.ylabel('psa', fontsize=18)
plt.show()

#%%
tab_usada['grupos_cluster_knn'] = cluster.fit(x).labels_
tab_usada.drop(["num_ro5_violations", "ro3_pass", "num_lipinski_ro5_violations"], axis=1).sort_values(by='grupos_cluster_knn')

#%%
sns.pairplot(tab_usada[['mw_freebase', 'psa', 'grupos_cluster_knn']], hue='grupos_cluster_knn')
