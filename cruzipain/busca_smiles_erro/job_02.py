#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import sqlite3
import numpy as np
from calcula_tanimoto import *

# In[2]:


import sys
sys.path.append("../../module/")
from lipinski import *


# In[3]:

con = sqlite3.connect(f"../../data/propriedades.db")
dataframe = pd.read_sql('select * from dados', con)
dataframe


# In[4]:

smiles = ''

data_com_tanimoto = calcula_tanimoto(smiles, dataframe, 20)


