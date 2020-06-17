#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import norm
get_ipython().run_line_magic('matplotlib', 'inline')


# In[3]:


HEAVY=pd.read_csv('heavy-cdr3-dist.csv')
LIGHT=pd.read_csv('light-cdr3-dist.csv')

HEAVY.head()


# In[16]:


plt.figure(figsize=(30,25))
sns.set(font_scale=1.5)
sns.set_style("white")

sns.kdeplot(HEAVY['heavy v'], color='#0f4c75', label='heavy', linewidth=8, shade=True, legend=False)
sns.kdeplot(LIGHT['lightv'], color='#bbe1fa', label='light', linewidth=8, shade=True, legend=False)


# plt.xlim(0,40)
# plt.ylim(0,0.2)

plt.savefig('COV2-CDR3.png')


# In[9]:


plt.figure(figsize=(30,30))
sns.set(font_scale=1.5)
sns.set_style("white")

sns.kdeplot(HEAVY['heavy v'], color='#0f4c75', label='heavy', linewidth=8, shade=True)

# plt.xlim(0,40)
# plt.ylim(0,0.2)

plt.savefig('COV2-CDR3-HC.png')


# In[10]:


plt.figure(figsize=(30,30))
sns.set(font_scale=1.5)
sns.set_style("white")

sns.kdeplot(LIGHT['lightv'], color='#bbe1fa', label='light', linewidth=8, shade=True)

# plt.xlim(0,40)
# plt.ylim(0,0.2)

plt.savefig('COV2-CDR3-LC.png')


# In[ ]:




