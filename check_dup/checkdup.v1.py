#!/usr/bin/env python
#!-*- coding:utf-8 -*-

"""Used to check depth
depth = depth - dup
"""

import  os
import sys
import colorama
import re

colorama.init()

Usage = """\033[1;3;32m
  used to check depth(depth - duplicate)
  Usage: checkdup infile  depth

\033[0m"""
 
if len(sys.argv) < 2:
    print(Usage)
    print(exit('need input file and design depth   !!!'))


def add_color(text, color=colorama.Fore.BLUE):
    return '%s %s %s'% (color, text, colorama.Fore.RESET)

def getpos(name, title):
    ntitle = [x.lower() for x in title]
    if name.lower() in ntitle:
        pos = ntitle.index(name.lower())
    else:
        print("%s is not in %s" % (name.lower(), ntitle))
        
    return pos

    
def checkdup(infile, depth):
    
    pattern = re.compile(r'.* [(](.*?)%[)]', re.S)    

    namelist = []
    duplist = []
    ori_depthlist = []
    with open(infile, 'r') as indata:
        for line in indata:
            lline =line.strip().split('\t')
            if line.startswith('ReportID'):
                namepos = getpos('Sample:', lline)
                duppos = getpos('Duplicate:', lline)
                deppos = getpos('Average_sequencing_depth_on_target:', lline)
            else:
                name = lline[namepos]
                dup = re.findall(pattern, lline[duppos])[0]            ## 129829842 (99.68%)
          
                ori_depth = lline[deppos] 
                namelist.append(name)
                duplist.append(dup)
                ori_depthlist.append(ori_depth)
            
        bad_depth_name = []     ## record sam which not stasify the depth
        
        n =1    ## used to cal sam number
        for sam,dup,ori_depth in zip(namelist, duplist, ori_depthlist):
            final_depth = float(ori_depth.strip("%")) - float(dup.strip("%"))
            print(add_color("\t".join((str(n), sam, dup, ori_depth, str(final_depth)))))
            n += 1

            if final_depth < float(depth.strip("%")):
                bad_depth_name.append([sam, dup, ori_depth, final_depth])
             
        
        print( "\n the sample below , depth is lower than %s\n" %  depth)
        print( add_color('%s\t%s\t%s\t%s' % ('#name', '#Duplicate', '#depth', '#depth-Duplicate'), colorama.Fore.RED) )
        for i in bad_depth_name:
            print(colorama.Back.YELLOW)
            print(i)
        
        
if __name__ == "__main__":
    infile = sys.argv[1]
    depth = sys.argv[2]
    
    checkdup(infile, depth)
    
    
