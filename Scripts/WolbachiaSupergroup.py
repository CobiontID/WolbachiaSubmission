from __future__ import division
import argparse
import configparser
import os
import sys
from ete3 import Tree

parser = argparse.ArgumentParser()
parser.add_argument("-b", type=str, action='store', dest='bin', metavar='INPUT',help='define bin name')
parser.add_argument("-f", type=str, action='store', dest='fam', metavar='INPUT',help='define family')
parser.add_argument("-s", type=str, action='store', dest='host', metavar='INPUT',help='define host scientific species name')
parser.add_argument("-t", type=str, action='store', dest='tree', metavar='INPUT',help='define iqtree SSU')
parser.add_argument("-ta", type=str, action='store', dest='tax', metavar='INPUT',help='define SSU tax file')
parser.add_argument("-c", type=str, action='store', dest='ctg', metavar='INPUT',help='define novel SSU ctg name')
parser.add_argument("-o", type=str, action='store', dest='out', metavar='OUT',help='define output filename')
parser.add_argument("-sts", type=str, action='store', dest='sts', metavar='STS',help='define stsfile')
parser.add_argument("-i", type=str, action='store', dest='tolid', metavar='INPUT',help='define tolid')
parser.add_argument("-na", type=str, action='store', dest='namesfile', metavar='NAMES',help='NCBI names.dmp')
parser.add_argument("-no", type=str, action='store', dest='nodesfile', metavar='NODES',help='NCBI nodes.dmp')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')
results = parser.parse_args()

def readNames(names_tax_file):
    '''
    input:
    - name.dmp (NCBI Taxonomy)
    output:
    - dictionary of form {node: name}
    - dictionary of form {sci name: node}
    '''
    tax_names = {}
    tax_names_reverse= {}
    syn_names={}
    with open(names_tax_file, 'r') as nodes_tax:
        for line in nodes_tax:
            node = [field.strip() for field in line.split('|')]
            if 'scientific' in line:
                tax_names[node[1]] = node[0]
                tax_names_reverse[node[0]] = node[1]
            if 'synonym' in line:
                syn_names[node[1]] = node[0]
    return tax_names_reverse,tax_names,syn_names

def readNodes(nodes_tax_file):

    '''
    input:
    - nodes.dmp (NCBI Taxonomy)
    output:
    - dictionary of form {parent: node}
    - dictionary of form {node: type}
    '''

    tax_nodes = {}
    tax_types = {}
    with open(nodes_tax_file, 'r') as nodes_tax:
        for line in nodes_tax:
            node = [field.strip() for field in line.split('|')]         #make list of line
            tax_nodes[node[0]] = node[1]                                #couple node with parent
            tax_types[node[0]] = node[2]                                #couple node with rank
    return tax_nodes, tax_types

def getTaxParent(tax_nodes, tax_types, taxid, ranking):

    '''
    input:
    - dictionary of form {parent: node} (readNodes output)
    - dictionary of form {node: type} (readNodes output)
    - taxid
    output:
    - dictionary of form {tax_id: [descendants]}
    '''
    tax_parents = {}
    node = str(taxid)
    if node not in tax_types:                                       #check if node in nodes.dmp
        tax_parents[node] = None
        print('[Warning] Could not find {} in nodes.dmp while parsing taxonomy hierarchy\n'.format(node))

    else:
        parent = tax_nodes[node]                                    #get parent for current node
        tax_parents[node] = [parent]                                #add node to dictionary

        while parent != tax_nodes[parent] and tax_types[parent]!= ranking:    #stop when parent = parent(parent) (i.e. 1 = 1)
            parent = tax_nodes[parent]                              #get parent of parent
            tax_parents[node].append(parent)                        #add parent to node in dictionary

    return tax_parents

genus=""
family=""
l=open(results.fam,'r')
for line in l:
    line=line.strip()
    genus=line.split(',')[0]
    family=line.split(',')[1]
l.close()

taxparents,taxtypes=readNodes(results.nodesfile)
taxnames,namestax,syn_names=readNames(results.namesfile)
hosttax=0
if results.host.replace("_"," ") in namestax:
    hosttax=namestax[results.host.replace("_"," ")]
elif results.host.replace("_"," ") in syn_names:
    hosttax=syn_names[results.host.replace("_"," ")]
lineage=getTaxParent(taxparents,taxtypes,hosttax,'order')
ordername=taxnames[lineage[hosttax][-1]]
print(lineage[hosttax][-1])
print(ordername)

multiple_host=False
orders_present=[]

matches={}
matches['Wolbachia'] = ['Arthropoda', 'Nematoda']
matches['Spiroplasma'] = ['Arthropoda']
matches['Cardinium'] = ['Arthropoda']
matches['Rickettsiella'] = ['Arthropoda']
matches['Arsenophonus'] = ['Arthropoda']
#matches['Lariskella'] = ['Arthropoda']
matches['Mesenet'] = ['Arthropoda']
#matches['Hepanticola'] = ['Arthropoda']
#matches['Rickettsia'] = ['Arthropoda']

#matches['Tremblaya'] = ['Pseudococcidae']  
matches['Buchnera'] = ['Aphididae'] 
#matches['Fukatsuia'] = ['Aphididae']
#matches['Regiella'] = ['Aphididae']
#matches['Hamiltonella'] = ['Aphididae']
#matches['Serratia'] = ['Aphididae']
#matches['Sulcia'] = ['Auchenorrhyncha']    
#matches['Hodgkinia'] = ['Auchenorrhyncha']
#matches['Baumannia'] = ['Auchenorrhyncha']
matches['Zinderia'] = ['Auchenorrhyncha']
matches['Carsonella'] = ['Psylloidea']
matches['Portiera'] = ['Aleyrodidae']
#matches['Nardonella'] = ['Curculionidae']
#matches['Walczuchella'] = ['Monophlebidae']
#matches['Uzinura'] = ['Diaspididae']

l=open(results.tax,'r')
for line in l:
    line=line.strip()
    if any([x in line for x in matches[genus]]) and ordername not in line:
        lastpart=line.split(';')[-2]
        if lastpart in namestax:
            lastparttax=namestax[lastpart]
        elif lastpart in syn_names:
            lastparttax=syn_names[lastpart]
        if lastparttax:
            if taxtypes[lastparttax] == 'order' and lastpart != ordername:
                if lastpart not in orders_present:
                    orders_present.append(lastpart)
                    multiple_host=True
            else:
                line_lineage=getTaxParent(taxparents,taxtypes,lastparttax,'order')
                line_order=taxnames[line_lineage[lastparttax][-1]]
                print(line_order)
                if line_order != 'root' and line_order != ordername:
                    orders_present.append(line_order)
                    multiple_host=True
l.close()

print(orders_present)

biospecimen=""
k=open(results.sts,'r')
for line in k:
    line=line.strip()
    if results.tolid in line:
        biospecimen=line.split('\t')[7]
        if biospecimen == '':
            biospecimen=line.split('\t')[6]
k.close()

o=open(results.out,'w')
if genus in ['Spiroplasma','Cardinium','Rickettsiella','Arsenophonus']:
    if multiple_host == False:
        o.write(genus+' endosymbiont of '+results.host.replace("_"," ")+'\tNovel Species\t'+results.host.replace("_"," ")+'\tPRJEB40665\tNovel endosymbionts from dToL samples\t'+results.bin+'\n')
    else:
        o.write(genus+' sp.\tNovel Species\t'+results.host.replace("_"," ")+'\tPRJEB40665\tNovel endosymbionts from dToL samples\t'+results.bin+'\n')
elif genus in ['Lariskella','Mesenet','Tremblaya']:
    if multiple_host == False:
        o.write('Candidatus '+genus+' endosymbiont of '+results.host.replace("_"," ")+'\tNovel Species\t'+results.host.replace("_"," ")+'\tPRJEB40665\tNovel endosymbionts from dToL samples\t'+results.bin+'\n')
    else:
        o.write('Candidatus '+genus+' sp.\tNovel Species\t'+results.host.replace("_"," ")+'\tPRJEB40665\tNovel endosymbionts from dToL samples\t'+results.bin+'\n')
elif genus == 'Sulcia':
    if multiple_host == False:
        o.write('Candidatus Karelsulcia muelleri ('+results.host.replace("_"," ")+')\tNovel Species\t'+results.host.replace("_"," ")+'\tPRJEB40665\tNovel endosymbionts from dToL samples\t'+results.bin+'\n')
    else:
        o.write('Candidatus Karelsulcia  muelleri\tNovel Species\t'+results.host.replace("_"," ")+'\tPRJEB40665\tNovel endosymbionts from dToL samples\t'+results.bin+'\n')       
elif genus == 'Buchnera':
    if multiple_host == False:
        o.write('Buchnera aphidicola ('+results.host.replace("_"," ")+')\tNovel Species\t'+results.host.replace("_"," ")+'\tPRJEB40665\tNovel endosymbionts from dToL samples\t'+results.bin+'\n')
    else:
        o.write('Buchnera aphidicola\tNovel Species\t'+results.host.replace("_"," ")+'\tPRJEB40665\tNovel endosymbionts from dToL samples\t'+results.bin+'\n') 
elif genus == 'Portiera':
    if multiple_host == False:
        o.write('Portiera aleyrodidarum endosymbiont of '+results.host.replace("_"," ")+'\tNovel Species\t'+results.host.replace("_"," ")+'\tPRJEB40665\tNovel endosymbionts from dToL samples\t'+results.bin+'\n')
    else:
        o.write('Portiera aleyrodidarum\tNovel Species\t'+results.host.replace("_"," ")+'\tPRJEB40665\tNovel endosymbionts from dToL samples\t'+results.bin+'\n') 
elif genus == 'Carsonella':
    if multiple_host == False:
        o.write('Candidatus Carsonella ruddii ('+results.host.replace("_"," ")+')\tNovel Species\t'+results.host.replace("_"," ")+'\tPRJEB40665\tNovel endosymbionts from dToL samples\t'+results.bin+'\n')
    else:
        o.write('Candidatus Carsonella ruddii\tNovel Species\t'+results.host.replace("_"," ")+'\tPRJEB40665\tNovel endosymbionts from dToL samples\t'+results.bin+'\n')
elif genus == 'Wolbachia':
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
        while len(node) < 3:
            #print(node)
            node = node.up

        sequences=0
        for leaf in node:
            if 'ptg' in leaf.name or 'atg' in leaf.name:
                sequences=sequences+1
        
        if sequences == len(node):
            node = node.up

        for leaf in node:
            #print(leaf.name)
            if leaf.name != ctg and  '_' in leaf.name:
                supergroup=leaf.name.split('_')[-1]
                #print(supergroup)
                if supergroup not in supergroups:
                    supergroups.append(supergroup)
        if len(supergroups) == 1:
            ctgs[ctg]=supergroups[0]
            if multiple_host == True:
                print(ctg+':Wolbachia sp. (group '+supergroups[0]+')'+'\t'+results.bin)
            else:
                print(ctg+':Wolbachia endosymbiont (group '+supergroups[0]+') of '+results.host.replace("_"," ")+'\t'+results.bin)
        else:
            ctgs[ctg]='Unclear'
            if multiple_host == True:
                print(ctg+':Wolbachia sp. '+'\t'+results.bin)
            else:
                print(ctg+':Wolbachia endosymbiont of '+results.host.replace("_"," ")+'\t'+results.bin)

    res = len(list(set(list(ctgs.values())))) == 1

    if res == True:
        supergroup_name=list(ctgs.values())[0]
        if supergroup_name != "Unclear":
            if multiple_host == True:
                o.write('Wolbachia sp. (group '+supergroup_name+') '+'\tNovel Species\t'+results.host.replace("_"," ")+'\tPRJEB40665\tNovel endosymbionts from dToL samples\t'+results.bin+'\n')
            else:
                o.write('Wolbachia endosymbiont (group '+supergroup_name+') of '+results.host.replace("_"," ")+'\tNovel Species\t'+results.host.replace("_"," ")+'\tPRJEB40665\tNovel endosymbionts from dToL samples\t'+results.bin+'\n')
        else:
            if multiple_host == True:
                o.write('Wolbachia sp. '+'\tNovel Species\t'+results.host.replace("_"," ")+'\tPRJEB40665\tNovel endosymbionts from dToL samples\t'+results.bin+'\n')
            else:
                o.write('Wolbachia endosymbiont of '+results.host.replace("_"," ")+'\tNovel Species\t'+results.host.replace("_"," ")+'\tPRJEB40665\tNovel endosymbionts from dToL samples\t'+results.bin+'\n')
    else:
        if multiple_host == True:
            o.write('Wolbachia sp. '+'\tNovel Species\t'+results.host.replace("_"," ")+'\tPRJEB40665\tNovel endosymbionts from dToL samples\t'+results.bin+'\n')
        else:
            o.write('Wolbachia endosymbiont of '+results.host.replace("_"," ")+'\tNovel Species\t'+results.host.replace("_"," ")+'\tPRJEB40665\tNovel endosymbionts from dToL samples\t'+results.bin+'\n')
o.close()
