from __future__ import division
import argparse
import configparser
import os
import sys
import json

parser = argparse.ArgumentParser()
parser.add_argument("-i", type=str, action='store', dest='input', metavar='INPUT',help='define taxid non-redundant file')
parser.add_argument("-o", type=str, action='store', dest='output', metavar='OUTPUT',help='define taxid output file with comma-separated bins')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')
results = parser.parse_args()

seen=[]
f =open(results.input,'r')
for record in f:
    idname=record.strip().split('\t')[0]
    idname_conv=idname.replace(' ','%20').strip()
    idname_conv1=idname.replace(' ','_').replace('(','_').replace(')','_').strip()
    #print(idname)
    cmd='curl -X GET --header "Accept: application/json" "https://www.ebi.ac.uk/ena/taxonomy/rest/scientific-name/'+idname_conv+'" --output '+idname_conv1+'.json'
    os.system(cmd)
    if os.path.getsize(idname_conv1+'.json') > 11:
        f = open(idname_conv1+'.json')
        data = json.load(f)
        if 'taxId' in data[0]:
            seen.append(idname)
            print(data[0]['taxId'])
        else:
            print('NO RECORD YET FOR '+idname)
        #os.system('rm sample.json')
    else:
        print('NO RECORD YET FOR '+idname)
    cmd='rm '+idname_conv1+'.json'
    os.system(cmd)
f.close()

f =open(results.input,'r')
o=open(results.output,'a')
for record in f:
    record=record.strip()
    idname=record.strip().split('\t')[0]
    if idname not in seen:
        o.write(record+'\n')
o.close()
f.close()