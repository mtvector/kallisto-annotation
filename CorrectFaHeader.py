####USAGE:
####python CorrectFaHeader.py GTF.gtf Input.fa Input.CorrectedForKallisto.fa -I
####Changes header and outputs proper tr2g.txt files for kallisto bus


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


GTF=dataframe(sys.argv[1])
GTFtrans=GTF.loc[GTF['feature'].isin(['transcript']),:]
#GTFtrans.index=GTFtrans['transcript_id']
GTFtrans.index=[re.sub('\.[0-9]+','',tid) for tid in GTFtrans['transcript_id']]
GTFtrans = GTFtrans.fillna('')
if len(sys.argv)>4:
    appendix=sys.argv[4]
else:
    appendix=''

counter = 0
newt2g=[[],[],[]]
with open(sys.argv[3],'w') as newfasta:
    with open(sys.argv[2],'r') as fasta:
        for line in fasta:
            if line[0] == '>':
                counter+=1
                tid=re.split('\s|:|\|',re.sub('>|\n','',line))[0]
                #To remove transcript version
                tid=re.sub('\.[0-9]+','',tid)
                newfasta.write('>'+tid+'.'+ str(counter)+ appendix+' gene_id:'+ GTFtrans.loc[tid,'gene_id']  + ' gene_name:' +GTFtrans.loc[tid,'gene_name'] +' '+ str(counter)+'\n')
                newt2g[0].append(tid+'.'+ str(counter)+ appendix)
                newt2g[1].append(GTFtrans.loc[tid,'gene_id'])
                newt2g[2].append(GTFtrans.loc[tid,'gene_name'])
            else:
                newfasta.write(line)

pd.DataFrame(newt2g).T.to_csv(os.path.join(os.path.dirname(sys.argv[2]),re.sub('\.fa','_t2g.txt',sys.argv[2])),sep='\t',header=False,index=False)
