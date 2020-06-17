#!/usr/bin/env python
# coding: utf-8

# In[4]:


import os
import csv


# In[1]:


count = {}
count_l={}

with open('cdr3s.csv') as fin, open('cdr3s-counted.csv', 'w') as fout:
    for line in fin:
        ls = line.strip().split(',')
        if ls[1] not in count:
            count[ls[1]] = 0
        count[ls[1]] += 1
        
    fin.seek(0,0)

    for line in fin:
        ls = line.strip().split(',')
        if ls[2] not in count_l:
            count_l[ls[2]] = 0
        count_l[ls[2]] += 1
        
    fin.seek(0,0)
    FIRST = True
    
    for line in fin:
        if FIRST:
            fout.write('mongo_id,heavy v,lightv,count-heavy,count-light\n')
            FIRST = False
        else:
            ls = line.strip().split(',')
            fout.write(line.strip() + ',' + str(count[ls[1]]) + ',' + str(count_l[ls[2]]) + '\n')


# In[ ]:




