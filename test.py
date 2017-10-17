
txt = open('list.txt').read()
import re
ls = re.split('[\n:\t]',txt)
lsnt = [e for e in ls if re.findall('.gz$',e) and re.findall('^nt',e)]
fo = open('list.txt','w')
for e in lsnt:
    fo.write('wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/'+e+'\n')
fo.close()