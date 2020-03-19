import pybedtools
from pybedtools import BedTool
import pandas as pd
import os
import sys
ifile=sys.argv[1]
#ifile='/wynton/group/ye/mtschmitz/refdata2/rhemac10/test_rhemac10/intron_sorted.bed'
efile=sys.argv[2]
#efile='/wynton/group/ye/mtschmitz/refdata2/rhemac10/test_rhemac10/transcript_sorted.bed'

i = BedTool(ifile)
e = BedTool(efile)
df=i.intersect(e, wa=False, wb=True).to_dataframe()
df.iloc[:,[0,1,2,6]].to_csv(os.path.join(os.path.dirname(ifile), 'intron_named.bed'),header=False,index=False,sep='\t')
