#!/usr/bin/env python
#-*- coding:utf-8 -*-

"""analysis hpa
"""

"""样本名称不能使用_分割，因为plink的利用vcf生成ped的过程中，Fid 和 SamID 的分割符是 _
"""

import os
import sys


#file with genename, one gene one row
def split_gene(genename):    
    code = """ ### split gene 
    for gene in `less  %s`; do
        echo ${gene} > ${gene}
    done """ % (genename)
    with open('step0.split_gene.sh', 'w') as outf:
        outf.write(code)
    # os.system(code)
    
def gene_list(genename):
    genelist = []
    with open(genename) as outdata:
        for line in outdata:
            line = line.strip()
            genelist.append(line)
        return genelist
    
    
##step1 : extract gene
def extract_gene(infile, genename):
    code = """ ##extract gene for infile
    for name in `less  %s `;do
        filter_gene -i  %s -o ${name}.filter.xls -gl ${name} 
    done""" % (genename, infile)
    with open('step1.extract_gene.sh', 'w') as outf:
        outf.write(code)
    os.system(code)
    
## get index
def getpos(name, title):
    ntitle = [x.lower() for x in title]
    if name.lower() in ntitle:
        pos = ntitle.index(name.lower())
        
    return pos


## get needed vcf file
def file2vcf(filterfile):
   
    with open(filterfile, 'r') as indata, open(filterfile.replace('xls', 'vcf'), 'w') as outf:
        for line in indata:
            title = line.strip().split('\t')
            if line.startswith('Priority'):
                #title = line.strip().split('\t')
                chr = getpos('CHROM', title)
               
                filter = getpos('FILTER', title)
                info = getpos('INFO', title)
                ori_ref = getpos('Ori_ref', title)
                outf.write('%s\t%s\n' % ('\t'.join(title[chr:filter+1]), '\t'.join(title[info, ori_ref]) ))
            elif title[3] == ".":
                title[2] = '_'.join((title[1], title[2], title[4], title[5]))
                outf.write('%s\t%s\n' % ('\t'.join(title[chr:filter+1]), '\t'.join(title[info, ori_ref]) ))
            else:
                outf.write('%s\t%s\n' % ('\t'.join(title[chr:filter+1]), '\t'.join(title[info, ori_ref]) ))
                
 
## plink analysis                
def plink(genename):
    code = """
    for name in `less  %s`;do
        plink --vcf ${name}.filter.vcf  --recode  --make-bed  --out ${out}.filter.vcf.plink
    done """ % genename
    with open('step2.plink.sh', 'w') as outf:
       outf.write(code)
       
       
## re_combine plink bed file

def plink2ped(plinkfile):
    with open(plinkfile, 'r') as indata, open(os.path.join(plinkfile, '.tmp'), 'w') as outdata:
        for line in indata:
            lline = line.strip().split('\t')
            string = '\t'.join(lline[6::]) + '\n'
            outdata.write(string)
          
            
## use sample_info to generate ped file

def info2ped(sample_info):
    with open(sample_info, 'r') as indata, open('sample_info_ped', 'w') as outdata:
        for line in indata:
            lline = line.strip('#').split('\t')
            if line.startswith('#') and 'sex' in [x.lower() for x in lline]:
                fam = getpos('familyid', lline)
                sam = getpos('sampleid', lline)
                sex = getpos('sex', lline)
                pn = getpos('normal/patient', lline)
                print(pn)
            elif not line.startswith('#'):
                if lline[sex] == "M":
                    Sex = 1
                elif lline[sex] == "F":
                    Sex = 2
                else:
                    Sex = 3
                
                if lline[pn] == "N":
                    PN = 1
                elif lline[pn] == "P":
                    PN = 2
                else:
                    PN = 0
                
                string = '\t'.join((lline[fam], lline[sam], '0', '0', str(Sex), str(PN)) ) + '\n'
                outdata.write(string)                                                
            else:
                print(line)
   
##  use *plink.bim to generate *HWE.info
def bim2info(genename):
    code = """
    for name in `less %s`;do
        paste sample_info_ped ${name}.filter.vcf.plink.ped.tmp > ${name}.filter.HWE.ped
        awk -F "\\t" -v OFS="\\t" '{if($2 != ".") print $2,$4; else print $4,$4}' > ${name}.filter.HWE.info
    done """ % genename
    with open('step3.addinfo.sh', 'w') as outf:
        outf.write(code)


## hpa analysis
def hpa(genename):
    code = """
    for name in `less %s`;do
        java -jar /NJPROJ1/DISEASE/WORK/zhaoyan/Haploview/Haploview.jar -pedfile ${name}.filter.HWE.ped  -info ${name}.filter.HWE.info \
        -png  -assocCC -maxDistance 5000    -blockoutput  GAB    -dprime -nogui
    done """ % genename
    with open('step4.hpa.sh', 'w') as outf:
       outf.write(code)
       
   
def main(genename, infile, sample_info):
    genelist = gene_list(genename)
    split_gene(genename)
    extract_gene(infile, genename)
    info2ped(sample_info)
    
    for gene in genelist:
        filterfile = os.path.join(gene, '.filter.xls')
        file2vcf(filterfile)
        
    plink(genename)
    
    for gene in genelist:
        plinkfile = os.path.join(gene, '.filter.vcf.plink.ped')
        plink2ped(plinkfile)
        
    bim2info(genename)
    
    hpa(genename)   
        
   
   
if __name__ == "__main__":
    genename = sys.argv[1]
    infile = sys.argv[2]
    sample_info = sys.argv[3]
    main(genename, infile, sample_info)