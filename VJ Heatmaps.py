#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
from pandas import DataFrame
import numpy as np
import os
import re
import statistics 
import scipy
from scipy.stats import zscore
import seaborn as sns
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')


# In[2]:


#the healthy and ebola heavies are flipped - but correct otherwise


# In[3]:


#read in all the gene distributions
HEAVY=pd.read_csv('heavy.csv')
KAPPA=pd.read_csv('kappa.csv')
LAMBDA=pd.read_csv('lambda.csv')

HEAVY.head(5)


# In[4]:


HEAVY_in=HEAVY.set_index('V')
KAPPA_in=KAPPA.set_index('V')
LAMBDA_in=LAMBDA.set_index('V')


HEAVY_in.head(5)


# In[7]:


fig, (ax1) = plt.subplots(1, 1, sharex=False, sharey=True, figsize = (160,30))
sns.set(font_scale=2.0)

sns.heatmap(HEAVY_in, vmin = -5, vmax= 5.0, cmap='coolwarm', robust=True, square=True, linewidths=2.0, cbar= True, cbar_kws={"fraction": 0.15, "shrink":1.0, "aspect":60}, ax=ax1)

fig.tight_layout()

plt.show()

fig.savefig('heatmapsHC.png')


# In[6]:


fig, (ax1, ax2) = plt.subplots(1, 2, sharex=False, sharey=False, figsize = (160,30))
sns.set(font_scale=2.0)

sns.heatmap(LAMBDA_in, vmin = -5, vmax= 5.0, cmap='coolwarm', robust=True, square=True, linewidths=2.0, cbar= False, ax=ax1)
sns.heatmap(KAPPA_in, vmin = -5, vmax= 5.0, cmap='coolwarm', robust=True, square=True, linewidths=2.0, cbar= True, cbar_kws={"fraction": 0.15, "shrink":1.0, "aspect":60}, ax=ax2)

fig.tight_layout()

plt.show()

fig.savefig('heatmapsLIGHTS.png')


# In[ ]:


# fig, (ax1, ax2) = plt.subplots(1, 2, sharex=False, sharey=True, figsize = (160,30))
# sns.set(font_scale=2.0)

# sns.heatmap(EBOV_lambda, vmin = -5, vmax= 5.0, cmap='coolwarm', robust=True, square=True, linewidths=2.0, cbar= False, ax=ax2)
# sns.heatmap(HEALTH_lambda, vmin = -5, vmax= 5.0, cmap='coolwarm', robust=True, square=True, linewidths=2.0, cbar= False, ax=ax5)

# fig.tight_layout()

# plt.show()

# fig.savefig('heatmapsKAPPA.png')

