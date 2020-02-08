#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import os
import warnings


# In[2]:


# 从筛选后图像保存路径中找到相应序号，回到表格中筛选
path_1 = r'D:\UCAS\adastro_once_used\QSO_redshift'
QSO_redshift = []
for root,dirs,files in os.walk(path_1):
    for file in files:
        QSO_redshift.append(int(os.path.splitext(file)[0]))
QSO_redshift.sort()

data_1 = pd.read_csv("flux_with_redshift.csv")
data_1['QSO'] = np.nan
for i in range(len(QSO_redshift)):
    warnings.filterwarnings('ignore')
    data_1['QSO'][QSO_redshift[i]] = 1
    
data_1.to_csv('flux_with_redshift_qso.csv', float_format="%.6f", index=0)


# In[3]:


path_2 = r'D:\UCAS\adastro_once_used\QSO_photoz'
QSO_photoz = []
for root,dirs,files in os.walk(path_2):
    for file in files:
        QSO_photoz.append(int(os.path.splitext(file)[0]))
QSO_photoz.sort()

data_2 = pd.read_csv("flux_with_photoz.csv")
data_2['QSO'] = np.nan
for i in range(len(QSO_photoz)):
    warnings.filterwarnings('ignore')
    data_2['QSO'][QSO_photoz[i]] = 1
    
data_2.to_csv('flux_with_photoz_qso.csv', float_format="%.6f", index=0)


# In[ ]:




