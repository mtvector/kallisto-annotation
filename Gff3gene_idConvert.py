import gffpandas
import gffpandas.gffpandas as gffpd
import sys
import re

annotation = gffpd.read_gff3(sys.argv[1])

newattr=[]
for x in annotation.df['attributes']:
    s=x.split('source_gene_common_name=')[1]
    s=s.split(';')[0]
    if s !='None':
        newattr.append(re.sub('gene_id=.+?(?=;)','gene_id='+s,x))
    else:
        newattr.append(x)    
        
annotation.df['attributes']=newattr
annotation.to_gff3(re.sub('.gff3','.geneid.gff3',sys.argv[1]))