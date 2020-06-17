#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
get_ipython().run_line_magic('matplotlib', 'inline')


# In[6]:


HEAVY=pd.read_csv('COV2-SHM-HC.csv')
HEAVY.head()


# In[20]:


g = sns.catplot(x='Chain', y='heavy_percent_id',
               data=HEAVY, kind="violin", linewidth=14, linecolor='black', palette=['#0f4c75'])
g.fig.set_figwidth(30)
g.fig.set_figheight(30)
g.set(ylim=(79, 105))

g.savefig("COV2SHM-HEAVY.png")


# In[13]:


LIGHT=pd.read_csv('COV2-SHM-LC.csv')
LIGHT.head()


# In[19]:


g = sns.catplot(x='Chain', y='heavy_percent_id',
               data=LIGHT, kind="violin", linewidth=14, linecolor='black', palette=['#bbe1fa'])
g.fig.set_figwidth(30)
g.fig.set_figheight(30)
g.set(ylim=(79, 105))

g.savefig("COV2SHM-LIGHT.png")

