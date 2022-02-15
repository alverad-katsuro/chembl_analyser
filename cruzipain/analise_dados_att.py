# %%
import pandas as pd
import sqlite3
import numpy as np

# %%
con = sqlite3.connect('../data/propriedades.db')
df = pd.read_sql('select * from dados', con)

# %%
cruz = pd.read_csv('data/cruzain_id_chemb.csv', sep=';')

# %%
df = df.drop(['cx_most_apka','cx_most_bpka','cx_logp','cx_logd'], axis=1)

# %%
df = df.drop(df.index[df.molecular_species.isnull()])

# %%
import seaborn as sns
from matplotlib import pyplot as plt
from sklearn.cluster import KMeans
from sklearn import preprocessing

# %%
pd.set_option('display.float_format', lambda x: '%.2f' % x)

# %%
df.query("canonical_smiles == 'c1c[nH]c(CN2CCc3ccccc3C2)c1'")

# %%
df.describe()

# %%
target = df.query("canonical_smiles == 'c1c[nH]c(CN2CCc3ccccc3C2)c1'").to_dict('record')

# %%
target = target[0]

# %%
target

# %%
axes = df.hist(bins=30,figsize=(15,15))
for i, ax in enumerate(axes.flatten()):
    title = ax.title.__dict__["_text"]
    if title:
        ax.axvline(x=target[title], color='r', linestyle='dashed', linewidth=2)

# %%
df.columns

# %%
rotulos = ['mw_freebase', 'alogp', 'psa', 'rtb', 'full_mwt', 'aromatic_rings', 'heavy_atoms', 'qed_weighted', 'mw_monoisotopic', 'hba_lipinski', 'hbd_lipinski']

# %%
def cria_query(propriedade_base: list, rotulos: dict, intervalo=0.3):
    string_query = ""
    for indice in range(len(rotulos)):
        maximo = propriedade_base[rotulos[indice]] * (1 + intervalo)
        minimo = propriedade_base[rotulos[indice]] * (1 - intervalo)
        if (indice == len(rotulos) - 1):
            string_query += f"{rotulos[indice]} < {maximo} and {rotulos[indice]} >= {minimo}"
        else:
            string_query += f"{rotulos[indice]} < {maximo} and {rotulos[indice]} >= {minimo} and "
    return string_query

# %%
df.query(cria_query(target, rotulos, 0.5))

# %%
#plt.hist(df.alogp, edgecolor='black',bins=30)
#plt.axvline(x=target[title], color='r', linestyle='dashed', linewidth=2)
#plt.title("AlogP")
#plt.grid()
#plt.show()

# %%
dir(df)

# %%
prep_clust = pd.get_dummies(df.query("num_lipinski_ro5_violations == 0"), columns=["molecular_species"])
prep_clust

# %%
prep_clust = prep_clust.drop(['chembl_id', "canonical_smiles"], axis=1)

# %%
df.query("chembl_id in @cruz['Molecule ChEMBL ID']").index

# %%
prep_clust

# %%
cruzain = df.query("chembl_id in @cruz['Molecule ChEMBL ID']")
cruzain

# %%
import sys
sys.path.append("../module/")
from calcula_tanimoto import *

# %%
cruzain_similaridade = []
jobs = cruzain.canonical_smiles.to_numpy().tolist()
for smile in jobs:
    cruzain_similaridade.append(calcula_tanimoto(smile, cruzain, 20))


# %%


# %%
clust_cru = pd.get_dummies(cruzain.drop(['chembl_id', "canonical_smiles"], axis=1), columns=["molecular_species"])
clust_cru

# %%


# %%
ab = preprocessing.normalize(clust_cru, norm='l1')
cluster_ab = KMeans(n_clusters=30, random_state=5)
grupos_ab = cluster_ab.fit(ab).labels_

# %%
plt.figure(figsize=(10,7))
plt.scatter(x=clust_cru.mw_freebase, y=clust_cru.psa, c=grupos_ab, s=50)
plt.xlabel('mw_freebase', fontsize=18)
plt.ylabel('psa', fontsize=18)
plt.show()

# %%
x = preprocessing.normalize(prep_clust, norm='l1')
cluster = KMeans(n_clusters=30, random_state=29)
cluster.fit(x).labels_

# %%
plt.figure(figsize=(10,7))
plt.scatter(x=prep_clust[['mw_freebase']], y=prep_clust.psa, c=cluster.fit(x).labels_, s=50)
plt.xlabel('mw_freebase', fontsize=18)
plt.ylabel('psa', fontsize=18)
plt.show()

# %%
plt.figure(figsize=(10,7))
plt.scatter(x=prep_clust.mw_freebase, y=prep_clust.hba_lipinski, c=cluster.fit(x).labels_, s=50)
plt.xlabel('alogp', fontsize=18)
plt.ylabel('psa', fontsize=18)
plt.show()

# %%
plt.figure(figsize=(10,7))
plt.scatter(x=prep_clust.alogp, y=prep_clust.psa, c=cluster.fit(x).labels_, s=50)
plt.xlabel('alogp', fontsize=18)
plt.ylabel('psa', fontsize=18)
plt.show()

# %%
plt.figure(figsize=(10,7))
ax = plt.axes(projection='3d')
# Data for a three-dimensional line
#zline = np.linspace(0, 15, 1000)
#xline = np.sin(zline)
#yline = np.cos(zline)
#ax.plot3D(xline, yline, zline, 'gray')

# Data for three-dimensional scattered points
ax.scatter3D(prep_clust.mw_freebase, prep_clust.psa, prep_clust.qed_weighted, c=cluster.fit(x).labels_);
plt.show

# %%
plt.figure(figsize=(10,7))
grp = plt.axes(projection='3d')
# Data for a three-dimensional line
#zline = np.linspace(0, 15, 1000)
#xline = np.sin(zline)
#yline = np.cos(zline)
#ax.plot3D(xline, yline, zline, 'gray')

# Data for three-dimensional scattered points
grp.scatter3D(prep_clust.mw_freebase, prep_clust.psa, prep_clust.alogp, c=cluster.fit(x).labels_);
plt.show

# %%
sns.pairplot(prep_clust[['mw_freebase', 'psa', 'grupos_cluster_knn']], hue='grupos_cluster_knn')


# %%
tab_usada['grupos_cluster_knn'] = cluster.fit(x).labels_
tab_usada.drop(["num_ro5_violations", "ro3_pass", "num_lipinski_ro5_violations"], axis=1).sort_values(by='grupos_cluster_knn')


