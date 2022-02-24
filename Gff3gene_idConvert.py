import gffpandas
import gffpandas.gffpandas as gffpd
import sys
import re
import bioinfokit
from bioinfokit.analys import gff
annotation = gffpd.read_gff3(sys.argv[1])
print('starting')
newattr=[]
for x in annotation.df['attributes']:
    s=x.split('Name=')[1]
    #s=x.split('source_gene_common_name=')[1]
    s=s.split(';')[0]
    print(s)
    if s !='None':
        newattr.append(re.sub('gene_id=.+?(?=;)','gene_id='+s,x))
    else:
        try:
            g=x.split('ID=')[1]
            g=g.split(';')[0]
            newattr.append(re.sub('=None','='+g,x))
            print(re.sub('=None','='+g,x))
        except:
            newattr.append(x)
        
annotation.df['attributes']=newattr
annotation.to_gff3(re.sub('.gff3','.geneid.gff3',sys.argv[1]))
gff.gff_to_gtf(file=re.sub('.gff3','.geneid.gff3',sys.argv[1]))
