#! /bin/python3
# For combining the "expected count" from individual result output into one file for DEseq analysis.
import sys
import pandas as pd
import re
l=sys.argv[1:]
df=pd.read_csv(l[0], sep='\t')
new=pd.DataFrame(df[['transcript_id','gene_id','expected_count']])
pattern=r'(.*)(\.rsem\.isoforms\.results)'
result=re.match(pattern,l[0])
#print(result.group(1))
new.rename(columns={'expected_count':result.group(1)},inplace=True)
#print(new)

for file in l[1:]:
    df=pd.read_csv(file, sep='\t')
    result=re.match(pattern,file)
    new[result.group(1)]=df['expected_count']
new.to_csv('isoform_count.txt',sep='\t',index=False)
