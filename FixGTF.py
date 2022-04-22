#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import sys
from collections import defaultdict
import gzip
import pandas as pd
import re
import csv

GTF_HEADER  = ['seqname', 'source', 'feature', 'start', 'end', 'score',
               'strand', 'frame']
R_SEMICOLON = re.compile(r'\s*;\s*')
R_COMMA     = re.compile(r'\s*,\s*')
R_KEYVALUE  = re.compile(r'(\s+|\s*=\s*)')


def dataframe(filename):
    """Open an optionally gzipped GTF file and return a pandas.DataFrame.
    """
    # Each column is a list stored as a value in this dict.
    result = defaultdict(list)

    for i, line in enumerate(lines(filename)):
        for key in line.keys():
            # This key has not been seen yet, so set it to None for all
            # previous lines.
            if key not in result:
                result[key] = [None] * i

        # Ensure this row has some value for each column.
        for key in result.keys():
            result[key].append(line.get(key, None))

    return pd.DataFrame(result)


def lines(filename):
    """Open an optionally gzipped GTF file and generate a dict for each line.
    """
    fn_open = gzip.open if filename.endswith('.gz') else open

    with fn_open(filename) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            else:
                yield parse(line)


def parse(line):
    """Parse a single GTF line and return a dict.
    """
    result = {}

    fields = line.rstrip().split('\t')

    for i, col in enumerate(GTF_HEADER):
        result[col] = _get_value(fields[i])

    # INFO field consists of "key1=value;key2=value;...".
    infos = [x for x in re.split(R_SEMICOLON, fields[8]) if x.strip()]

    for i, info in enumerate(infos, 1):
        # It should be key="value".
        try:
            key, _, value = re.split(R_KEYVALUE, info, 1)
        # But sometimes it is just "value".
        except ValueError:
            key = 'INFO{}'.format(i)
            value = info
        # Ignore the field if there is no value.
        if value:
            result[key] = _get_value(value)

    return result


def _get_value(value):
    if not value:
        return None

    # Strip double and single quotes.
    value = value.strip('"\'')

    # Return a list if the value has a comma.
    if ',' in value:
        value = re.split(R_COMMA, value)
    # These values are equivalent to None.
    elif value in ['', '.', 'NA']:
        return None

    return value


# In[10]:


common_name='Name'
common_id='gene_id'
transcript_or_gene='transcript'
gtf = dataframe('/wynton/group/ye/mtschmitz/refdata2/rhemac10/CAT_chang/Rhesus.gtf')
gtf = gtf.fillna('')
print(gtf)
print(list(gtf.columns))


# In[12]:


tgtf=gtf.loc[gtf['feature']=='transcript',:]
tgtf.loc[:,['transcript_id',common_id]]
tgtf.drop_duplicates(subset = ['transcript_id'], keep = 'first', inplace = True) 

ggtf=gtf.loc[gtf['feature']=='transcript',:]
ggtf.drop_duplicates(subset = [common_id], keep = 'first', inplace = True) 
gdict=dict(zip(ggtf[common_id],ggtf[common_name]))
#ggtf.loc[:,[common_id,common_name]]
tgtf[common_name]=list(tgtf[common_id].replace(gdict))


# In[16]:


tgtf[common_name]=list(tgtf[common_id].replace(gdict))


# In[13]:


gtf['feature'].unique()


# In[14]:


ggtf


# In[15]:


tgtf


# In[ ]:




