#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os,sys,argparse,re

Usage = """
\033[1;3;32m
计算2个软件的交集作为新的gff结果

overlap ：  根据模式，自己指定
\033[0m"""

if len(sys.argv) < 4:
    print("python gff_intersect.py infile1.gff infile2.gff outgff overlap")
    print(Usage)
    exit()



 
def safe_open(filename, mode='r'):
    try:
        if mode == 'w':
            dirname = os.path.dirname(filename)
            if dirname and not os.path.exists(dirname):
                os.makedirs(dirname)
        if filename.endswith('.gz'):
            mode += 'b'
            import gzip
            return gzip.open(filename, mode)
        return open(filename, mode)
    except Exception as e:
        exit(e)

def get_pos(name,title,name2=''):
    ntitle = [i.lower() for i in title]
    if name.lower() in ntitle:
        pos = ntitle.index(name.lower())
    elif name2 != '' and name2.lower() in ntitle:
        pos = ntitle.index(name2.lower())
    else:
        if name2 != '':
           exit('Wrong: %s and %s not in title.' %(name,name2))
        else:
           exit('Wrong: %s not in title.' %name)
    return pos

def gff2bed(gff,bed):
    infile=safe_open(gff,"r")
    outfile=safe_open(bed,"w")
    for everyline in infile:
        everyline = everyline.strip().split('\t')
        if 'Chr' in everyline:
            print(everyline)
            chrpos = get_pos('chr', everyline)
            startpos = get_pos('start', everyline)
            endpos = get_pos('end', everyline)
            typepos = get_pos('svtype', everyline)
            #qFuncpos = get_pos('Func.refGene', everyline)
        else:
            chr1  = everyline[chrpos]
            start=everyline[startpos]
            end=everyline[endpos]
            svtype = everyline[typepos]   #突变类型
            #Func = everyline[Funcpos]     ##作为featurexinx
            #detail = everyline[8]
            #anno = "TCHR=na;TSTART=na;SVID=" + str(SVID)
            if svtype in  ["DEL","DUP"]:
                outfile.write("\t".join([chr1,start,end])+"\n")
        

# def write_bed(infile):
#     outfile = infile + '.bed'
#     return gff2bed(infile, outfile)


def intersect(bed1, bed2, cc):
    ##通过intersectBed计算bed区间的重叠区域
    intesct_tmp = 'share_intersect_%s_sv.temp' % cc
    code = '''# find instersct of two software : > %f length
    intersectBed -a %s -b %s -f %f -r  -wa -wb > %s
    '''% (float(cc),bed1,bed2,float(cc),intesct_tmp)
    return os.system(code)


def file2dict(infile):
    """将输入文件，转化为字典，
    key为chr-start-end
    value为整行内容
    """
    tmpdict = {}
    with open(infile, 'r') as indata:
        for line in indata:
            lline = line.strip().split('\t')
            if 'Chr' in lline:
                chrpos = get_pos('chr', lline)
                startpos = get_pos('start', lline)
                endpos = get_pos('end', lline)
            else:
                chr = lline[chrpos]
                start = lline[startpos]
                end = lline[endpos]
                tmpdict['-'.join([chr, start, end])] = line
                
    return tmpdict
                

def findkey(keys, dicts, out):
    """查找key， 如果key存在，输出到文件out中
    """
    with open(out, 'w') as outf:
        for i in keys:
            if i in dicts:
                #print(dicts[i])
                outf.write(dicts[i])
            
            

def intect_file(infile1, infile2, cc):
    """根据输出的tmp文件，输出在输入文件中的内容
    """
    tmp = 'share_intersect_%s_sv.temp' % cc
    key1 = []  ##infile1 的查找key
    key2 = []  ##infile2 的查找key
    
    with open(tmp, 'r') as indata:
        for line in indata:
            lline = line.strip().split('\t')
            key1.append('-'.join(lline[:3]))
            key2.append('-'.join(lline[3:]))
    dict1 = file2dict(infile1)
    dcit2 = file2dict(infile2)
    
    out1 = infile1.split('.')[0] + '_sect_with_' + cc + '.xls'
    out2 = infile2.split('.')[0] + '_sect_with_' + cc + '.xls'
    
    findkey(key1, dict1, out1)
    findkey(key2, dcit2, out2)
    
def main():
    
    ##file 2 bed
    bed1 = infile1 + '.bed'
    bed2 = infile2 + '.bed'
    
    gff2bed(infile1, bed1)
    gff2bed(infile2, bed2)
     
    ## get intersect bed, out is tmp
    intersect(bed1, bed2, cc)

    ## file 2 dict
    dict1 = file2dict(infile1)
    dict2 = file2dict(infile2)
    
    ## write alone result
    intect_file(infile1, infile2, cc)
    

if __name__ == "__main__":
    infile1 = sys.argv[1]  #软件1的gff文件
    infile2 = sys.argv[2]  #软件2的gff文件
    suffix = sys.argv[3]   #交集结果
    cc = sys.argv[4]       #overlap的数值
    main()
    
    





















# ##按照ID对文件进行排序
# with safe_open(intesct_tmp,"r") as denovotemp, safe_open(intesct_file,"w") as denovoFi:
#                 templist=denovotemp.readlines()
#                 sortlist=sorted(templist,key=lambda item:int(re.search("SVID=(.*?)\t",item).group(1)))
#                 for line in sortlist:
#                     denovoFi.write(line)

# ###注释
# type = 'SV'
# Mtype=type+"Type"
# MID=type+"ID"


# if type=="SV":
#     cmd='''
# awk -F"\\t" -v OFS="\\t" '{print $1,$6,$4,$2,$3,".",".",".","Size="$3-$2+1";"$5";%s="$4}' %s > %s
# ''' %(Mtype,intesct_file,intesct_gff)
# # if ref=="hg19":
# #     cmd+='''
# # sh /PUBLIC/software/HUMAN/ANNOVAR_2017Jul16/Var_annotation_disease_ANNOVAR2017Jul16_v4.6.sh -t %s %s %s
# # '''% (Mtype,denovogff,childID)
# # elif ref=="hg38":
# #     cmd+='''
# # sh /PUBLIC/software/HUMAN/ANNOVAR_2015Dec14/Var_annotation_disease_ANNOVAR2015Dec14_hg38.sh -t %s %s %s
# # '''% (Mtype,denovogff,childID)
# # cmd+='''
# # rm -f %s %s %s
# # '''% (familydir+"/*.temp",familydir+"/*.itx",familydir+"/*.ctx")

# os.system(cmd)

