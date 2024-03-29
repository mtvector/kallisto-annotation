import os
import sys
from collections import defaultdict
import gzip
import pandas as pd
import re
import csv

common_name='gene_name'
gene_id='gene_id'
transcript_or_gene='gene'

if len(sys.argv)>4:
    common_name=sys.argv[4]#'Name'
if len(sys.argv)>3:
    gene_id=sys.argv[3]#'gene_id'
if len(sys.argv)>2:
    transcript_or_gene=sys.argv[2]#'transcript'

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

gtf = dataframe(sys.argv[1])
gtf = gtf.fillna('')
print(gtf)
print(list(gtf.columns))
tgtf=gtf.loc[gtf['feature']=='transcript',:]
tgtf.loc[:,['transcript_id',gene_id]]
tgtf.drop_duplicates(subset = ['transcript_id'], keep = 'first', inplace = True) 

print(gtf)
#gtf=gtf.loc[gtf['feature']=='exon',:]
column_subset=gtf.columns[8:]
tid=list(gtf['transcript_id'])
print(gtf.loc[:,column_subset])
gtf.loc[:,column_subset] = '"' + gtf.loc[:,column_subset].astype(str) + '"'
gtf.loc[:,column_subset] = gtf.loc[:,column_subset].astype(str) + ';'
for header in column_subset:
    gtf[header] = header + ' ' + gtf[header]
gtf=gtf.fillna('')
gtf['features'] = gtf[gtf.columns[8:]].apply(lambda x: ' '.join(x), axis=1)
gtf['transcript_id']=tid

bed=gtf.loc[:,['seqname','start','end','transcript_id','score','strand','source','score','features']]
bed['start']=bed['start'].astype(int)
bed.loc[:,'start']=bed.loc[:,'start']-1
bed.to_csv(re.sub('gtf','bed',sys.argv[1]), sep='\t', header=False, index=False, quoting=csv.QUOTE_NONE)

print(gtf['feature'])
print('boolean',gtf['feature']==transcript_or_gene)
ggtf=gtf.loc[gtf['feature']==transcript_or_gene,:]
ggtf.drop_duplicates(subset = [gene_id], keep = 'first', inplace = True) 
print('ggtf',ggtf)
gdict=dict(zip(ggtf[gene_id],ggtf[common_name]))
#ggtf.loc[:,[gene_id,common_name]]
tgtf[common_name]=list(tgtf[gene_id].replace(gdict))
tgtf.loc[:,['transcript_id',gene_id,common_name]].to_csv(os.path.join(os.path.dirname(sys.argv[1]),'tr2g.txt'),header=False,index=False,sep='\t')

