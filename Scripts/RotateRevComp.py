from __future__ import division
import argparse
import configparser
import os
from symtable import symtable
import sys

parser = argparse.ArgumentParser()
parser.add_argument("-gff", type=str, action='store', dest='gff',help='define gff file')
parser.add_argument("-fa", type=str, action='store', dest='fa',help='define fasta file')
parser.add_argument("-o", type=str, action='store', dest='outfile',help='define output file name')
parser.add_argument("-f", type=str, action='store', dest='fam',help='define family')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')
args = parser.parse_args()

species=args.outfile.split('/')[-1].split('.fa')[0]
fastafile = args.fa

if args.fam == 'Anaplasmataceae':
    genename='hemE'
elif args.fam in ['Flavobacteriales','Betaproteobacteria','Erwiniaceae','Halomonadaceae']:
    genename='dnaE'
else:
    genename='dnaA'

#double-check dnaE for Tremblaya, Buchnera, Sulcia, Carsonella, Portiera
#Fukatsuia, Uzinura and Walczuchella

seqname=""
seqlen=0
k=open(fastafile,'r')
for line in k:
    line=line.strip()
    if line.startswith('>'):
        seqname=line.split('>')[1]
    else:
        seqlen+=int(len(line))
print(fastafile)

l=open(args.gff,'r')
for record in l:
    record=record.strip()
    if record.startswith(seqname):
        start=record.split('\t')[3]
        stop=record.split('\t')[4]
        ctg=record.split('\t')[0]
        strand=record.split('\t')[6]
        if genename in record or genename+'_1;' in record:
            print(start+'\t'+stop+'\t'+strand)
            if strand == '+':
                n=open(species+'.HemE.end.bed','w')
                n.write(seqname+'\t'+'0'+'\t'+str(start)+'\n')
                n.close()
                n=open(species+'.HemE.start.bed','w')
                n.write(seqname+'\t'+str(start)+'\t'+str(int(seqlen)+1)+'\n')
                n.close()
                cmd='seqtk subseq '+fastafile +' '+species+".HemE.start.bed | grep -v '>' | tr -d '\n' > "+species+'.HemE.start.fa'
                os.system(cmd)
                cmd='seqtk subseq '+fastafile +' '+species+".HemE.end.bed | grep -v '>' | tr -d '\n' > "+species+'.HemE.end.fa'
                os.system(cmd)
                cmd='cat '+species+'.HemE.start.fa '+species+'.HemE.end.fa | fold > '+species+'.HemE.fa'
                os.system(cmd)
                cmd='(echo ">'+seqname+'" && cat '+species+'.HemE.fa) > '+args.outfile
                print(cmd)
                os.system(cmd)
                cmd='rm '+species+'.HemE.start* '+species+'.HemE.end* '+species+'.HemE.fa'
                os.system(cmd)
            elif strand == '-':
                n=open(species+'.HemE.start.bed','w')
                n.write(seqname+'\t'+'0'+'\t'+str(int(stop)-1)+'\n')
                n.close()
                n=open(species+'.HemE.end.bed','w')
                n.write(seqname+'\t'+str(stop)+'\t'+str(int(seqlen)+1)+'\n')
                n.close()
                cmd='seqtk subseq '+fastafile +' '+species+".HemE.start.bed | seqtk seq -r | grep -v '>' | tr -d '\n' > "+species+'.HemE.start.fa'
                os.system(cmd)
                cmd='seqtk subseq '+fastafile +' '+species+".HemE.end.bed | seqtk seq -r | grep -v '>' | tr -d '\n' > "+species+'.HemE.end.fa'
                os.system(cmd)
                cmd='cat '+species+'.HemE.start.fa '+species+'.HemE.end.fa | fold > '+species+'.HemE.fa'
                os.system(cmd)
                cmd='(echo ">'+seqname+'" && cat '+species+'.HemE.fa) > '+args.outfile
                print(cmd)
                os.system(cmd)
                cmd='rm '+species+'.HemE.start* '+species+'.HemE.end* '+species+'.HemE.fa'
                os.system(cmd)                            

l.close()