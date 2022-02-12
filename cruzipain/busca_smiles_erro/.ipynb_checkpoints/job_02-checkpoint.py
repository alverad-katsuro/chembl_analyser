#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import sqlite3
import numpy as np


# In[2]:


import sys
sys.path.append("../../module/")
from lipinski import *


# In[3]:


from rdkit import Chem
# To add a new cell, type '# %%'
# To add a new markdown cell, type '# %% [markdown]'
# %%
import pandas as pd
import sqlite3
import seaborn as sns
from matplotlib import pyplot as plt
from sklearn.cluster import KMeans
from sklearn import preprocessing

from matplotlib import pyplot as ptl
import seaborn as sns
import sys, os
from rdkit import Chem
from rdkit.Chem import Crippen
from rdkit.Chem import Lipinski
from rdkit.Chem import Descriptors
from rdkit.Chem.QED import properties

import subprocess
import threading
import queue


# In[3]:


con = sqlite3.connect(f"../../data/propriedades.db")
dataframe = pd.read_sql('select * from dados', con)
dataframe


# In[4]:


teste = dataframe['canonical_smiles']


# In[5]:


x = teste.to_numpy().tolist()


# In[6]:


def ver_lip(x, i):
    file = open(f"smiles_erro_log_{i}.log",'a')
    for smiles in x:
        try:
            verifica_lipinski(smiles)
        except:
            file.write(f"{smiles}\n")
    file.close()


# In[7]:


def verifica_kuku(array, threads_num):
  jobs = np.array_split(array, threads_num)
  processos = []    
  for i in range(1, threads_num):
    processThread = threading.Thread(target=ver_lip, args=(jobs.pop(), i))
    processThread.start()
    processos.append(processThread)
  ver_lip(jobs.pop(), 0)
  for proc in processos:
    proc.join()


# In[ ]:


verifica_kuku(x,20)


# In[ ]:




