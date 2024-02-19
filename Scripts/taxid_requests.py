from __future__ import division
import argparse
import configparser
import os
import sys

parser = argparse.ArgumentParser()
parser.add_argument("-i", type=str, action='store', dest='input', metavar='INPUT',help='define taxid non-redundant file')
parser.add_argument("-o", type=str, action='store', dest='output', metavar='OUTPUT',help='define taxid output file with comma-separated bins')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')
results = parser.parse_args()

seen={}
f =open(results.input,'r')
for record in f:
    record=record.strip()
    string=str(';'.join(record.split('\t')[0:4]))
    if string not in seen:
        seen[string]=[]
    seen[string].append(record.split('\t')[5])

o=open(results.output,'a')
for taxs in seen:
    o.write(taxs.replace(';','\t')+'\t'+','.join(seen[taxs])+'\n')
o.close()