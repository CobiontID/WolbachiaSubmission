from __future__ import division
import argparse
import configparser
import os
import sys
from ete3 import Tree

parser = argparse.ArgumentParser()
parser.add_argument("-b", type=str, action='store', dest='bin', metavar='INPUT',help='define bin name')
parser.add_argument("-s", type=str, action='store', dest='host', metavar='INPUT',help='define host scientific species name')
parser.add_argument("-t", type=str, action='store', dest='tree', metavar='INPUT',help='define iqtree SSU')
parser.add_argument("-c", type=str, action='store', dest='ctg', metavar='INPUT',help='define novel SSU ctg name')
parser.add_argument("-o", type=str, action='store', dest='out', metavar='OUT',help='define output filename')
parser.add_argument("-sts", type=str, action='store', dest='sts', metavar='STS',help='define stsfile')
parser.add_argument("-i", type=str, action='store', dest='tolid', metavar='INPUT',help='define tolid')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')
results = parser.parse_args()

biospecimen=""
k=open(results.sts,'r')
for line in k:
    line=line.strip()
    if results.tolid in line:
        biospecimen=line.split('\t')[7]
        if biospecimen == '':
            biospecimen=line.split('\t')[6]
k.close()

ctgs={}
m =open(results.ctg,'r')
for record in m:
    if record.startswith('>'):
        ctg=record.split()[0].split('>')[1]
        ctgs[ctg]=""

t = Tree(results.tree)
for ctg in ctgs:
    supergroups=[]
    node = t.search_nodes(name=ctg)[0]
    i=0
    while i < 3:
        node = node.up
        i=i+1
    for leaf in node:
        if leaf.name != ctg:
            supergroup=leaf.name.split('_')[-1]
            if supergroup not in supergroups:
                supergroups.append(supergroup)
    if len(supergroups) == 1:
        ctgs[ctg]=supergroups[0]
        print(ctg+':Wolbachia endosymbiont (group '+supergroups[0]+') of '+results.host.replace("_"," ")+'\t'+results.bin)
    else:
        ctgs[ctg]='Unclear'
        print(ctg+':Wolbachia endosymbiont of '+results.host.replace("_"," ")+'\t'+results.bin)

res = len(list(set(list(ctgs.values())))) == 1

o=open(results.out,'w')
if res == True:
    supergroup_name=list(ctgs.values())[0]
    if supergroup_name != "Unclear":
        o.write('Wolbachia endosymbiont (group '+supergroup_name+') of '+results.host.replace("_"," ")+'\t'+results.bin+"\t"+biospecimen+"\n")
    else:
        o.write('Wolbachia endosymbiont of '+results.host.replace("_"," ")+'\t'+results.bin+'\t'+biospecimen+"\n")
else:
    o.write('Wolbachia endosymbiont of '+results.host.replace("_"," ")+'\t'+results.bin+"\t"+biospecimen+"\n")