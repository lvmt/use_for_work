#!/usr/bin/env python
#!-*- coding:utf-8 -*-

"""Used to check depth
depth = depth - dup
"""

import  os
import sys
sys.path.append('/ifs/TJPROJ3/DISEASE/share/Software/python/site-packages')
import colorama
import re
import argparse
import sys

colorama.init()

Usage = """\033[1;3;32m
  used to check depth(depth - duplicate) 
  py infile --depth  100

  used to check ori_depth
  py infile --dup 100 -oo 
 
  used to check datasize(10G or other)
  py infile --datasize 15
  
  used to check dup rate (20 or other)
  python infile --dup 20
  
  Usage: checkdup infile  --xx xx

\033[0m"""

if len(sys.argv) < 2:
    print(Usage)
    print(exit('\033[1;3;34m need input file and check Type   !!!'))

parser = argparse.ArgumentParser(description="check threshold", prog="checkfile", usage="%(prog)s [option]")

parser.add_argument("infile", help="inputfile")
parser.add_argument("--depth", help = "check the Threshold depth; the depth  (depth - dup)")
parser.add_argument('--datasize', help = "check the Threshold datasize")
parser.add_argument('--dup',  help = "check the Threshold depth")
parser.add_argument('--20x', help= "check the Threshold 20X ")
parser.add_argument("--oo", help="if just check the ori depth; store_true", action='store_true' )

argv =parser.parse_args()


## add color
def add_color(text, color=colorama.Fore.BLUE):
    return '%s %s %s'% (color, text, colorama.Fore.RESET)

## get pos
def getpos(name, title):
    ntitle = [x.lower() for x in title]
    if name.lower() in ntitle:
        pos = ntitle.index(name.lower())
    else:
        print("%s is not in %s" % (name.lower(), ntitle))
    return pos


## define all check list

namelist = []      
duplist = []
ori_depthlist = []
fra_20xlist = []
Maplist = []
pro_maplist = []
bedlist = []          ## 记录捕获效率
coveragelist = []     # 覆盖度

## check depth
pattern = re.compile(r'.* [(](.*?)%[)]', re.S)

def checkfile(infile):

    with open(infile, 'r') as indata:
        for line in indata:
            lline =line.strip().split('\t')
            if line.startswith('ReportID'):
                namepos = getpos('Sample:', lline)
                duppos = getpos('Duplicate:', lline)
                deppos = getpos('Average_sequencing_depth_on_target:', lline)
                fra_20xpos = getpos('Fraction_of_target_covered_with_at_least_20x:', lline)
                mappos = getpos('Mapped:', lline)
                pro_mappos = getpos('Properly mapped:', lline)
                bedpos = getpos('Fraction_of_effective_bases_on_target:', lline)
                coveragepos = getpos('Coverage_of_target_region:', lline)
                
            else:
                name = lline[namepos]
                dup = re.findall(pattern, lline[duppos])[0]            ## 129829842 (99.68%)
                ori_depth = lline[deppos]
                fra_20x = lline[fra_20xpos]
                mapp = lline[mappos]
                pro_map = lline[pro_mappos]
                bed = lline[bedpos]
                cov = lline[coveragepos]
                
                namelist.append(name)
                duplist.append(dup)
                ori_depthlist.append(ori_depth)
                fra_20xlist.append(fra_20x) 
                Maplist.append(mapp)
                pro_maplist.append(pro_map)
                bedlist.append(bed)
                coveragelist.append(cov)
          
      
## check depth
def checkdepth(depth):
    
    bad_depth_name = []     ## record sam which not stasify the depth(depth - dup
    bad_ori_depth_name = []   ## record sam which not stasify depth

    # n =1    ## used to cal sam number
    for sam,dup,ori_depth in zip(namelist, duplist, ori_depthlist):
        final_depth = float(ori_depth.strip("%")) - float(dup.strip("%"))

        if final_depth < float(depth.strip("%")):
            bad_depth_name.append([sam, dup, ori_depth, final_depth])

        if float(ori_depth.strip("%")) < float(depth):
            bad_ori_depth_name.append([sam, dup, ori_depth])

    alllist = zip(namelist, duplist, ori_depthlist)
        
    return(bad_depth_name, bad_ori_depth_name, alllist)


## check datasize
def checksize(infile, datasize):
    namelist = []
    rawdatalist = []
    with open(infile, 'r') as indata:
        for line in indata:
            lline = line.strip().split('\t')
            if line.startswith('Sample'):
                sampos = getpos('Sample name', lline)
                rawpos = getpos('Raw data(G)', lline)
            else:
                sam = lline[sampos]
                rawdata = lline[rawpos]
                namelist.append(sam)
                rawdatalist.append(rawdata)

    bad_rawdata_name = []

    n = 1

    for name, rawdata in zip(namelist, rawdatalist):
        print(add_color("\t".join((str(n), name, rawdata))))
        if float(rawdata) < float(datasize):
            bad_rawdata_name.append([name, rawdata])

        n += 1

    print('\n the sample below, rawdatasize is less than %s G\n' % datasize)
    print(add_color('%s\t%s\t%s\t' % ('#ID' ,'#Name', '#Datasize'), colorama.Fore.RED))

    m = 1
    for i in bad_rawdata_name:
        print(colorama.Back.GREEN)
        print(m, i)
        m += 1

## check dup ratw
def checkdup(dup):
    """判断dup是否超过一定的比率：
    输出的内容： ID，sam， depth， dup
    """
    bad_dup_list = []
    allset = zip(namelist, ori_depthlist, duplist)
    for sam, depth, tdup in allset:
        tdup = re.findall(pattern,dup)[0]            ## 129829842 (99.68%)
        if float(tdup) > float(dup):
            bad_dup_list.append([sam, depth, tdup])
    
    ## print result 
    n = 1
    for name, ori_depth, dup in zip(namelist, ori_depthlist, duplist):
        print(add_color('\t'.join((str(n), name, ori_depth, dup))))
        n += 1
    
    ## print bad dup sam
    print('\n the sample below, Dup is more than %s %\n' % dup)
    print(add_color('%s\t%s\t%s\t%s' % ('#ID' ,'#Name', '#Depth', '#Dup'), colorama.Fore.RED))
    
    m = 1
    for i in bad_dup_list:
        print(colorama.Back.GREEN)
        print(m, i)
        

def main():
    infile = argv.infile
    if argv.depth:
        checkfile(infile)
        depth = argv.depth
        bad_depth_name,bad_ori_depth_name, alllist = checkdepth(depth)
        
        if argv.oo:
            n = 1
            for sam,dup,ori_depth in alllist:
                print(add_color('\t'.join((str(n), sam, dup, ori_depth)) ))
                n +=1 
                
            m = 1
            print('\nthe sample below , depth is lower than %sX' % depth)
            print(add_color("\n%s\t%s\t%s\t%s" % ('#ID', '#Name', '#Dup', '#Ori_Depth'), colorama.Fore.RED))
            print(colorama.Back.GREEN)
            for i in bad_ori_depth_name:
                print(m, i)
                m += 1
                
        elif not argv.oo:
            k = 1
            for sam,dup,ori_depth in alllist:
                final_depth = float(ori_depth.strip("%")) - float(dup.strip("%"))
                print(add_color('\t'.join((str(k), sam, dup, ori_depth, str(final_depth) )) ))
                k +=1 
                
            j = 1
            print('\nthe sample below , (depth - dup) is lower than %sX' % depth)
            print(add_color("\n%s\t%s\t%s\t%s\t%s" % ('#ID', '#Name', '#Dup', '#Ori_depth', '#Depth - Dup'), colorama.Fore.RED))
            print(colorama.Back.GREEN)
            for i in bad_depth_name:
                print(j, i)
                j += 1
    
    if argv.dup:
        infile = argv.infile
        checkfile(infile)
        checkdup(dup)
        
    
        
    
    if argv.datasize:
        datasize = argv.datasize
        checksize(infile, datasize)

if __name__ == "__main__":
    main()
    

