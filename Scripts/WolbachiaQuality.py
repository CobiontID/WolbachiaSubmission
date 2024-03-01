from __future__ import division
import argparse
import configparser
import os
import sys
from datetime import datetime

parser = argparse.ArgumentParser()
parser.add_argument("-d", type=str, action='store', dest='dir', metavar='INPUT',help='define dir')
parser.add_argument("-d2", type=str, action='store', dest='dir2', metavar='INPUT',help='define orig dir')
parser.add_argument("-o", type=str, action='store', dest='out', metavar='INPUT',help='define output dir')
parser.add_argument("-t", type=str, action='store', dest='tolid', metavar='INPUT',help='define tolid')
parser.add_argument("-l", type=str, action='store', dest='binlist', metavar='INPUT',help='define binlist')
parser.add_argument("-g", type=str, action='store', dest='gfa', metavar='GFA',help='define gfa file')
parser.add_argument("-f", type=str, action='store', dest='fam', metavar='FAM',help='define family')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')
results = parser.parse_args()

genus=""
family=""
l=open(results.fam,'r')
for line in l:
    line=line.strip()
    genus=line.split(',')[0]
    family=line.split(',')[1]
l.close()

requirements={}
requirements['Anaplasmataceae']={}
requirements['Anaplasmataceae']['buscoscore']=96
requirements['Anaplasmataceae']['minlen']=700000
requirements['Anaplasmataceae']['maxlen']=2400000

requirements['Spiroplasmataceae']={}
requirements['Spiroplasmataceae']['buscoscore']=60
requirements['Spiroplasmataceae']['minlen']=700000
requirements['Spiroplasmataceae']['maxlen']=2400000

requirements['Amoebophilaceae']={}
requirements['Amoebophilaceae']['buscoscore']=35
requirements['Amoebophilaceae']['minlen']=700000
requirements['Amoebophilaceae']['maxlen']=1500000

requirements['Rickettsiaceae']={}
requirements['Rickettsiaceae']['buscoscore']=96
requirements['Rickettsiaceae']['minlen']=1200000
requirements['Rickettsiaceae']['maxlen']=2800000

requirements['Coxiellaceae']={}
requirements['Coxiellaceae']['buscoscore']=80
requirements['Coxiellaceae']['minlen']=1200000
requirements['Coxiellaceae']['maxlen']=2000000

requirements['Morganellaceae']={}
requirements['Morganellaceae']['buscoscore']=80
requirements['Morganellaceae']['minlen']=700000
requirements['Morganellaceae']['maxlen']=5000000

requirements['Midichloriaceae']={}
requirements['Midichloriaceae']['buscoscore']=90
requirements['Midichloriaceae']['minlen']=700000
requirements['Midichloriaceae']['maxlen']=2000000

requirements['Flavobacteriales']={}
requirements['Flavobacteriales']['buscoscore']=15
requirements['Flavobacteriales']['minlen']=100000
requirements['Flavobacteriales']['maxlen']=500000

requirements['Betaproteobacteria']={}
requirements['Betaproteobacteria']['buscoscore']=2
requirements['Betaproteobacteria']['minlen']=100000
requirements['Betaproteobacteria']['maxlen']=300000

requirements['Oxalobacteraceae']={}
requirements['Oxalobacteraceae']['buscoscore']=20
requirements['Oxalobacteraceae']['minlen']=100000
requirements['Oxalobacteraceae']['maxlen']=500000

requirements['Erwiniaceae']={}
requirements['Erwiniaceae']['buscoscore']=50
requirements['Erwiniaceae']['minlen']=300000
requirements['Erwiniaceae']['maxlen']=800000

requirements['Halomonadaceae']={}
requirements['Halomonadaceae']['buscoscore']=5
requirements['Halomonadaceae']['minlen']=200000
requirements['Halomonadaceae']['maxlen']=500000

binlistfile=open(results.binlist,'w')

m =open(results.dir,'r')
totalfound=0
total_busco=0
complete=0
dirshort=results.dir.split('/buscoAssembly/completeness_per_contig.txt')[0]
for rec in m:
    if not rec.startswith('#'):
        ctg=rec.split('\t')[0]
        totalfound+=int(rec.split('\t')[1])
        total_busco=int(rec.split('\t')[2])
        if float(rec.split('\t')[3].split('%')[0]) > requirements[family]['buscoscore'] and float(rec.split('\t')[5]) > requirements[family]['minlen'] and float(rec.split('\t')[5]) < requirements[family]['maxlen']:
            complete+=1
m.close()

totalfound2=0
total_busco2=0
complete2=0
ctgs_busco=[]
k=open(dirshort+'/busco/completeness_per_contig.txt')
for rec in k:
    if not rec.startswith('#'):
        ctg=rec.split('\t')[0]
        ctgs_busco.append(ctg)
        totalfound2+=int(rec.split('\t')[1])
        total_busco2=int(rec.split('\t')[2])
        if float(rec.split('\t')[3].split('%')[0]) > requirements[family]['buscoscore'] and float(rec.split('\t')[5]) > requirements[family]['minlen'] and float(rec.split('\t')[5]) < requirements[family]['maxlen']:
            complete2+=1
k.close()
print(str(totalfound2)+','+str(total_busco2)+','+str(complete2)+','+str(complete))

contignumber=1
m =open(results.dir,'r')
if (totalfound/total_busco >0.97 and round(totalfound/total_busco)==complete):
    for record in m:
        record=record.strip()
        if not record.startswith('#'):
            if float(record.split('\t')[3].split('%')[0]) > requirements[family]['buscoscore'] and float(record.split('\t')[5]) > requirements[family]['minlen'] and float(record.split('\t')[5]) < requirements[family]['maxlen']:
                print(record)
                ctg=record.split('\t')[0]
                if ctg.endswith('c'):
                    print(ctg)
                    gfafile=dirshort+'/hifiasm/hifiasm.p_ctg.gfa'
                    n=open(gfafile,'r')
                    coverage=0
                    for line in n:
                        line=line.strip()
                        if line.startswith('S\t'+ctg):
                            coverage=int(line.split('\t')[4].split(':')[2])
                    n.close()
                    #print(coverage)
                    today=datetime.today().strftime('%Y%m%d')
                    outdir=results.out+'/'+results.tolid+'.'+genus+'_sp_'+str(contignumber)+'.1.'+today
                    cmd="mkdir "+outdir
                    os.system(cmd)
                    y=open(outdir+'/'+results.tolid+'.'+genus+'_sp_'+str(contignumber)+'.1.yaml', 'w')
                    y.write('---'+'\n')
                    #y.write('species: '+cobiontname+'\n')
                    y.write('specimen: '+results.tolid+'.'+genus+'_sp_'+str(contignumber)+'\n')
                    y.write('projects:'+'\n')
                    y.write('  - darwin'+'\n')
                    y.write('data_location: Sanger RW'+'\n')
                    y.write('chromosome_list: '+outdir+'/'+results.tolid+'.'+genus+'_sp_'+str(contignumber)+'.1.tsv'+'\n')
                    y.write('cobiont_status: cobiont'+'\n')
                    y.write('primary: '+outdir+'/'+results.tolid+'.'+genus+'_sp_'+str(contignumber)+'.1.fa.gz'+'\n')
                    y.write('jira_queue: DS'+'\n')
                    #y.write('biosample: '+biosample+'\n')
                    y.write('coverage: '+str(coverage)+'\n')
                    #y.write('taxid: '+str(taxid)+'\n')
                    y.write('assembly_source: MarkerScan'+'\n')
                    y.write('pipeline:'+'\n')
                    y.write('  - MarkerScan (v1.0)'+'\n')
                    y.write('  - hifiasm (v0.14)'+'\n')
                    y.close()
                    lst=open(outdir+'/'+results.tolid+'.'+genus+'_sp_'+str(contignumber)+'.1.tsv', 'w')
                    lst.write(ctg+'\t1\tCircular-Chromosome\n')
                    lst.close()
                    lst2=open(dirshort+'/'+results.tolid+'.'+genus+'_sp_'+str(contignumber)+'.1.list', 'w')
                    lst2.write(ctg+'\t1\tCircular-Chromosome\n')
                    lst2.close()
                    fafile=open(outdir+'/'+results.tolid+'.'+genus+'_sp_'+str(contignumber)+'.1.fa', 'w')
                    fafile2=open(dirshort+'/'+results.tolid+'.'+genus+'_sp_'+str(contignumber)+'.1.fa', 'w')
                    fasta=dirshort+'/hifiasm/hifiasm.p_ctg.fasta'
                    n=open(fasta,'r')
                    inline=False
                    for line in n:
                        if line.startswith('>'):
                            if ctg in line:
                                inline=True
                                fafile.write(line)
                                fafile2.write(line)
                            else:
                                inline=False
                        else:
                            if inline == True:
                                fafile.write(line)
                                fafile2.write(line)
                    n.close()
                    fafile.close()
                    fafile2.close()
                    cmd="gzip "+outdir+'/'+results.tolid+'.'+genus+'_sp_'+str(contignumber)+'.1.fa'
                    os.system(cmd)
                    binlistfile.write(results.tolid+'.'+genus+'_sp_'+str(contignumber)+'.1\n')
                    contignumber+=1
                else:
                    #print(ctg)
                    ctgs=[]
                    ctgs.append(ctg)
                    #buscotable=dirshort+'/buscoAssembly/busco/run_rickettsiales_odb10/full_table.tsv'
                    buscotable=dirshort+'/buscoAssembly/full_table.tsv'
                    buscofound=[]
                    buscofound_additional={}
                    listfile=False
                    b=open(buscotable)
                    for line in b:
                        line=line.strip()
                        if ctg in line:
                            buscofound.append(line.split('\t')[0])
                        elif not line.startswith('#') and not 'Missing' in line:
                            #print(line)
                            if not line.split('\t')[2] in buscofound_additional:
                                buscofound_additional[line.split('\t')[2]]=[]
                            buscofound_additional[line.split('\t')[2]].append(line.split('\t')[0]) 
                    b.close()

                    b=open(buscotable)
                    for novelctg in buscofound_additional:
                        totalbusco=0
                        for line in b:
                            line=line.strip()
                            if novelctg in line:
                                if not line.split('\t')[0] in buscofound:
                                    #print(line)
                                    totalbusco+=1
                        print(novelctg+'\t'+str(totalbusco)+'\t'+str(len(buscofound_additional[novelctg]))+'\t'+','.join(buscofound_additional[novelctg]))
                        if totalbusco == len(buscofound_additional[novelctg]):
                            print(','.join(buscofound_additional[novelctg]))
                            ctgs.append(novelctg)
                    b.close()

                    today=datetime.today().strftime('%Y%m%d')
                    outdir=results.out+'/'+results.tolid+'.'+genus+'_sp_'+str(contignumber)+'.1.'+today
                    cmd="mkdir "+outdir
                    os.system(cmd)

                    gfafile=dirshort+'/hifiasm/hifiasm.p_ctg.gfa'
                    n=open(gfafile,'r')
                    final_coverage=0
                    if len(ctgs) == 1:
                        for line in n:
                            line=line.strip()
                            if line.startswith('S\t'+ctg):
                                final_coverage=int(line.split('\t')[4].split(':')[2])
                        print(final_coverage)
                        listfile=True
                        lst=open(outdir+'/'+results.tolid+'.'+genus+'_sp_'+str(contignumber)+'.1.list', 'w')
                        lst.write(ctg+'\t1\tLinear-Chromosome\n')
                        lst.close()
                        lst2=open(dirshort+'/'+results.tolid+'.'+genus+'_sp_'+str(contignumber)+'.1.list', 'w')
                        lst2.write(ctg+'\t1\tLinear-Chromosome\n')
                        lst2.close()
                    else:
                        calc_cov=0
                        total_len=0
                        for line in n:
                            line=line.strip()
                            if line.startswith('S') and line.split('\t')[1] in ctgs:
                                cov=int(line.split('\t')[4].split(':')[2])
                                if cov == 0:
                                    cov=1
                                length=int(line.split('\t')[3].split(':')[2])
                                total_len+=length
                                calc_cov+=(cov*length)
                                print(line.split('\t')[1]+'\t'+str(length)+'\t'+str(cov))    
                        final_coverage=int(calc_cov/total_len)
                    n.close()

                    y=open(outdir+'/'+results.tolid+'.'+genus+'_sp_'+str(contignumber)+'.1.yaml', 'w')
                    y.write('---'+'\n')
                    #y.write('species: '+cobiontname+'\n')
                    y.write('specimen: '+results.tolid+'.'+genus+'_sp_'+str(contignumber)+'\n')
                    y.write('projects:'+'\n')
                    y.write('  - darwin'+'\n')
                    y.write('data_location: Sanger RW'+'\n')
                    if listfile == True:
                        y.write('chromosome_list: '+outdir+'/'+results.tolid+'.'+genus+'_sp_'+str(contignumber)+'.1.tsv'+'\n')
                    y.write('cobiont_status: cobiont'+'\n')
                    y.write('primary: '+outdir+'/'+results.tolid+'.'+genus+'_sp_'+str(contignumber)+'.1.fa.gz'+'\n')
                    y.write('jira_queue: DS'+'\n')
                    #y.write('biosample: '+biosample+'\n')
                    y.write('coverage: '+str(final_coverage)+'\n')
                    #y.write('taxid: '+str(taxid)+'\n')
                    y.write('assembly_source: MarkerScan'+'\n')
                    y.write('pipeline:'+'\n')
                    y.write('  - MarkerScan (v1.0)'+'\n')
                    y.write('  - hifiasm (v0.14)'+'\n')
                    y.close()

                    fafile=open(outdir+'/'+results.tolid+'.'+genus+'_sp_'+str(contignumber)+'.1.fa', 'w')
                    fafile2=open(dirshort+'/'+results.tolid+'.'+genus+'_sp_'+str(contignumber)+'.1.fa', 'w')
                    fasta=dirshort+'/hifiasm/hifiasm.p_ctg.fasta'
                    n=open(fasta,'r')
                    inline=False
                    for line in n:
                        if line.startswith('>'):
                            if line.split()[0].split('>')[1] in ctgs:
                                print(line)
                                inline=True
                                fafile.write(line)
                                fafile2.write(line)
                            else:
                                inline=False
                        else:
                            if inline == True:
                                fafile.write(line)
                                fafile2.write(line)
                    n.close()
                    fafile.close()
                    fafile2.close()
                    cmd="gzip "+outdir+'/'+results.tolid+'.'+genus+'_sp_'+str(contignumber)+'.1.fa'
                    os.system(cmd)
                    binlistfile.write(results.tolid+'.'+genus+'_sp_'+str(contignumber)+'.1\n')
                    contignumber+=1
elif (totalfound2/total_busco2 >0.97 and round(totalfound2/total_busco2)==complete2):
    k=open(dirshort+'/busco/completeness_per_contig.txt')
    for record in k:
        record=record.strip()
        if not record.startswith('#'):
            if float(record.split('\t')[3].split('%')[0]) > 97 and float(record.split('\t')[5]) > 700000 and float(record.split('\t')[5]) < 2400000:
                print(record)
                ctg=record.split('\t')[0]
                if ctg.endswith('c'):
                    print(ctg)
                    gfafile=results.gfa
                    if results.gfa.endswith('gz'):
                        if not os.path.isfile(dirshort+'/'+results.gfa.split('/')[-1].split('.gz')[0]):
                            cmd='cp '+results.gfa+' '+dirshort+'/'
                            os.system(cmd)
                            #print(cmd)
                            cmd='gunzip '+dirshort+'/'+results.gfa.split('/')[-1]
                            #print(cmd)
                            os.system(cmd)
                        gfafile=dirshort+'/'+results.gfa.split('/')[-1].split('.gz')[0]
                    n=open(gfafile,'r')
                    coverage=0
                    for line in n:
                        line=line.strip()
                        if line.startswith('S\t'+ctg):
                            coverage=int(line.split('\t')[4].split(':')[2])
                    n.close()
                    #print(coverage)
                    today=datetime.today().strftime('%Y%m%d')
                    outdir=results.out+'/'+results.tolid+'.'+genus+'_sp_'+str(contignumber)+'.1.'+today
                    cmd="mkdir "+outdir
                    os.system(cmd)
                    y=open(outdir+'/'+results.tolid+'.'+genus+'_sp_'+str(contignumber)+'.1.yaml', 'w')
                    y.write('---'+'\n')
                    #y.write('species: '+cobiontname+'\n')
                    y.write('specimen: '+results.tolid+'.'+genus+'_sp_'+str(contignumber)+'\n')
                    y.write('projects:'+'\n')
                    y.write('  - darwin'+'\n')
                    y.write('data_location: Sanger RW'+'\n')
                    y.write('chromosome_list: '+outdir+'/'+results.tolid+'.'+genus+'_sp_'+str(contignumber)+'.1.tsv'+'\n')
                    y.write('cobiont_status: cobiont'+'\n')
                    y.write('primary: '+outdir+'/'+results.tolid+'.'+genus+'_sp_'+str(contignumber)+'.1.fa.gz'+'\n')
                    y.write('jira_queue: DS'+'\n')
                    #y.write('biosample: '+biosample+'\n')
                    y.write('coverage: '+str(coverage)+'\n')
                    #y.write('taxid: '+str(taxid)+'\n')
                    y.write('assembly_source: MarkerScan'+'\n')
                    y.write('pipeline:'+'\n')
                    y.write('  - MarkerScan (v1.0)'+'\n')
                    y.write('  - hifiasm (v0.14)'+'\n')
                    y.close()
                    lst=open(outdir+'/'+results.tolid+'.'+genus+'_sp_'+str(contignumber)+'.1.tsv', 'w')
                    lst.write(ctg+'\t1\tCircular-Chromosome\n')
                    lst.close()
                    lst2=open(dirshort+'/'+results.tolid+'.'+genus+'_sp_'+str(contignumber)+'.1.list', 'w')
                    lst2.write(ctg+'\t1\tCircular-Chromosome\n')
                    lst2.close()
                    fafile=open(outdir+'/'+results.tolid+'.'+genus+'_sp_'+str(contignumber)+'.1.fa', 'w')
                    fafile2=open(dirshort+'/'+results.tolid+'.'+genus+'_sp_'+str(contignumber)+'.1.fa', 'w')
                    #fasta=dirshort+'/hifiasm/hifiasm.p_ctg.fasta'
                    cmd='cp '+results.gfa.split('.noseq.gfa')[0]+'.fa.gz '+dirshort+'/'
                    os.system(cmd)
                    #print(cmd)
                    fasta=dirshort+'/'+results.gfa.split('/')[-1].split('.noseq.gfa')[0]+'.fa'
                    if not os.path.isfile(fasta):
                        cmd='gunzip '+dirshort+'/'+results.gfa.split('/')[-1].split('.noseq.gfa')[0]+'.fa.gz'
                        #print(cmd)
                        os.system(cmd)
                    n=open(fasta,'r')
                    inline=False
                    for line in n:
                        if line.startswith('>'):
                            if ctg in line:
                                inline=True
                                fafile.write(line)
                                fafile2.write(line)
                            else:
                                inline=False
                        else:
                            if inline == True:
                                fafile.write(line)
                                fafile2.write(line)
                    n.close()
                    fafile.close()
                    fafile2.close()
                    cmd="gzip "+outdir+'/'+results.tolid+'.'+genus+'_sp_'+str(contignumber)+'.1.fa'
                    os.system(cmd)
                    binlistfile.write(results.tolid+'.'+genus+'_sp_'+str(contignumber)+'.1\n')
                    contignumber+=1
                else:
                    #print(ctg)
                    ctgs=[]
                    ctgs.append(ctg)
                    #buscotable=dirshort+'/buscoAssembly/busco/run_rickettsiales_odb10/full_table.tsv'
                    buscotable=dirshort+'/busco/full_table.tsv'
                    buscofound=[]
                    buscofound_additional={}
                    listfile=False
                    b=open(buscotable)
                    for line in b:
                        line=line.strip()
                        if ctg in line:
                            buscofound.append(line.split('\t')[0])
                        elif not line.startswith('#') and not 'Missing' in line:
                            #print(line)
                            if not line.split('\t')[2] in buscofound_additional:
                                buscofound_additional[line.split('\t')[2]]=[]
                            buscofound_additional[line.split('\t')[2]].append(line.split('\t')[0]) 
                    b.close()

                    b=open(buscotable)
                    for novelctg in buscofound_additional:
                        totalbusco=0
                        for line in b:
                            line=line.strip()
                            if novelctg in line:
                                if not line.split('\t')[0] in buscofound:
                                    #print(line)
                                    totalbusco+=1
                        print(novelctg+'\t'+str(totalbusco)+'\t'+str(len(buscofound_additional[novelctg]))+'\t'+','.join(buscofound_additional[novelctg]))
                        if totalbusco == len(buscofound_additional[novelctg]):
                            print(','.join(buscofound_additional[novelctg]))
                            ctgs.append(novelctg)
                    b.close()

                    today=datetime.today().strftime('%Y%m%d')
                    outdir=results.out+'/'+results.tolid+'.'+genus+'_sp_'+str(contignumber)+'.1.'+today
                    cmd="mkdir "+outdir
                    os.system(cmd)

                    gfafile=dirshort+'/hifiasm/hifiasm.p_ctg.gfa'
                    n=open(gfafile,'r')
                    final_coverage=0

                    gfafile=results.gfa
                    if results.gfa.endswith('gz'):
                        if not os.path.isfile(dirshort+'/'+results.gfa.split('/')[-1].split('.gz')[0]):
                            cmd='cp '+results.gfa+' '+dirshort+'/'
                            os.system(cmd)
                            #print(cmd)
                            cmd='gunzip '+dirshort+'/'+results.gfa.split('/')[-1]
                            #print(cmd)
                            os.system(cmd)
                        gfafile=dirshort+'/'+results.gfa.split('/')[-1].split('.gz')[0]
                    n=open(gfafile,'r')

                    if len(ctgs) == 1:
                        for line in n:
                            line=line.strip()
                            if line.startswith('S\t'+ctg):
                                final_coverage=int(line.split('\t')[4].split(':')[2])
                        print(final_coverage)
                        listfile=True
                        lst=open(outdir+'/'+results.tolid+'.'+genus+'_sp_'+str(contignumber)+'.1.list', 'w')
                        lst.write(ctg+'\t1\tLinear-Chromosome\n')
                        lst.close()
                        lst2=open(dirshort+'/'+results.tolid+'.'+genus+'_sp_'+str(contignumber)+'.1.list', 'w')
                        lst2.write(ctg+'\t1\tLinear-Chromosome\n')
                        lst2.close()
                    else:
                        calc_cov=0
                        total_len=0
                        for line in n:
                            line=line.strip()
                            if line.startswith('S') and line.split('\t')[1] in ctgs:
                                cov=int(line.split('\t')[4].split(':')[2])
                                if cov == 0:
                                    cov=1
                                length=int(line.split('\t')[3].split(':')[2])
                                total_len+=length
                                calc_cov+=(cov*length)
                                print(line.split('\t')[1]+'\t'+str(length)+'\t'+str(cov))    
                        final_coverage=int(calc_cov/total_len)
                    n.close()

                    y=open(outdir+'/'+results.tolid+'.'+genus+'_sp_'+str(contignumber)+'.1.yaml', 'w')
                    y.write('---'+'\n')
                    #y.write('species: '+cobiontname+'\n')
                    y.write('specimen: '+results.tolid+'.'+genus+'_sp_'+str(contignumber)+'\n')
                    y.write('projects:'+'\n')
                    y.write('  - darwin'+'\n')
                    y.write('data_location: Sanger RW'+'\n')
                    if listfile == True:
                        y.write('chromosome_list: '+outdir+'/'+results.tolid+'.'+genus+'_sp_'+str(contignumber)+'.1.tsv'+'\n')
                    y.write('cobiont_status: cobiont'+'\n')
                    y.write('primary: '+outdir+'/'+results.tolid+'.'+genus+'_sp_'+str(contignumber)+'.1.fa.gz'+'\n')
                    y.write('jira_queue: DS'+'\n')
                    #y.write('biosample: '+biosample+'\n')
                    y.write('coverage: '+str(final_coverage)+'\n')
                    #y.write('taxid: '+str(taxid)+'\n')
                    y.write('assembly_source: MarkerScan'+'\n')
                    y.write('pipeline:'+'\n')
                    y.write('  - MarkerScan (v1.0)'+'\n')
                    y.write('  - hifiasm (v0.14)'+'\n')
                    y.close()

                    fafile=open(outdir+'/'+results.tolid+'.'+genus+'_sp_'+str(contignumber)+'.1.fa', 'w')
                    fafile2=open(dirshort+'/'+results.tolid+'.'+genus+'_sp_'+str(contignumber)+'.1.fa', 'w')
                    cmd='cp '+results.gfa.split('.noseq.gfa')[0]+'.fa.gz '+dirshort+'/'
                    os.system(cmd)
                    #print(cmd)
                    fasta=dirshort+'/'+results.gfa.split('/')[-1].split('.noseq.gfa')[0]+'.fa'
                    if not os.path.isfile(fasta):
                        cmd='gunzip '+dirshort+'/'+results.gfa.split('/')[-1].split('.noseq.gfa')[0]+'.fa.gz'
                        #print(cmd)
                        os.system(cmd)
                    n=open(fasta,'r')
                    inline=False
                    for line in n:
                        if line.startswith('>'):
                            if line.split()[0].split('>')[1] in ctgs:
                                print(line)
                                inline=True
                                fafile.write(line)
                                fafile2.write(line)
                            else:
                                inline=False
                        else:
                            if inline == True:
                                fafile.write(line)
                                fafile2.write(line)
                    n.close()
                    fafile.close()
                    fafile2.close()
                    cmd="gzip "+outdir+'/'+results.tolid+'.'+genus+'_sp_'+str(contignumber)+'.1.fa'
                    os.system(cmd)
                    binlistfile.write(results.tolid+'.'+genus+'_sp_'+str(contignumber)+'.1\n')
                    contignumber+=1               
else:
    binlistfile.write(results.tolid+'.'+genus+'_sp_'+str(contignumber)+'.1\n')
    fafile=dirshort+'/'+results.tolid+'.'+genus+'_sp_'+str(contignumber)+'.1.fa'
    orig_fafile=results.dir2+'/'+family+'/'+family+'.finalassembly.fa'
    if os.path.getsize(orig_fafile):
        cmd="cp "+orig_fafile+" "+fafile
        os.system(cmd)
    elif os.path.getsize(results.dir2+'/'+family+'/'+family+'.ctgs.fa'):
        k=open(results.dir2+'/'+family+'/'+family+'.ctgs.fa','r')
        o=open(fafile,'w')
        for line in k:
            line=line.strip
            if line.startswith('>'):
                if line.split('>')[1] in ctgs_busco:
                    found=1
                    o.write(line)
                else:
                    found=0
            else:
                if found == 1:
                    o.write(line)
        k.close()
        o.close()

    ctglist=[]
    n=open(orig_fafile,'r')
    for line in n:
        line=line.strip()
        if line.startswith('>'):
            ctglist.append(line.split('>')[1])
    n.close()

    calc_cov=0
    total_len=0
    if results.gfa.endswith('gz'):
        gfafile2=results.gfa.split('.p_ctg')[0]+'.a_ctg.noseq.gfa.gz'
        if not os.path.isfile(dirshort+'/'+results.gfa.split('/')[-1].split('.gz')[0]):
            cmd='cp '+results.gfa+' '+dirshort+'/'
            os.system(cmd)
            #print(cmd)
            cmd='gunzip '+dirshort+'/'+results.gfa.split('/')[-1]
            #print(cmd)
            os.system(cmd)
        cmd='cp '+gfafile2+' '+dirshort+'/'
        #print(cmd)
        os.system(cmd)
        cmd='gunzip '+dirshort+'/'+gfafile2.split('/')[-1]
        #print(cmd)
        os.system(cmd)
        gfafiles=[dirshort+'/'+results.gfa.split('/')[-1].split('.gz')[0], dirshort+'/'+gfafile2.split('/')[-1].split('.gz')[0]]
    else:
        gfafile2=results.gfa.split('.p_ctg')[0]+'.a_ctg.noseq.gfa'
        gfafiles=[results.gfa, gfafile2]
        
    for gfa in gfafiles:
        l=open(gfa,'r')
        for line in l:
            line=line.strip()
            if line.startswith('S') and line.split('\t')[1] in ctglist:
                cov=int(line.split('\t')[4].split(':')[2])
                if cov == 0:
                    cov=1
                length=int(line.split('\t')[3].split(':')[2])
                total_len+=length
                calc_cov+=(cov*length)
                print(line.split('\t')[1]+'\t'+str(length)+'\t'+str(cov))
        l.close()
    final_coverage=0
    if total_len > 0:
        final_coverage=int(calc_cov/total_len)
        manifestfile=open(dirshort+'/'+results.tolid+'.'+genus+'_sp_'+str(contignumber)+'.1.manifest.txt','w')
        manifestfile.write('STUDY\t\nSAMPLE\t\nASSEMBLYNAME\t'+results.tolid+'.'+genus+'_sp_'+str(contignumber)+'.1\nASSEMBLY_TYPE\tMetagenome-Assembled Genome (MAG)\nCOVERAGE\t'+str(final_coverage)+'X\nPROGRAM\tHifiasm\nPLATFORM\tPacBio Sequel II (HiFi)\nMOLECULETYPE\tgenomic DNA\nFASTA\t'+results.tolid+'.'+genus+'_sp_'+str(contignumber)+'.1.mag.fa\n')
        manifestfile.close()

    contignumber+=1
binlistfile.close()
m.close()