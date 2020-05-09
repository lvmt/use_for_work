#!/usr/bin/env python
#-*- coding:utf-8 -*-

import re 
import sys
import textwrap

def cp(ossfile, outfile, outdir):

    with open(ossfile, 'r') as fr , open(outfile, 'w') as fo:
        for line in fr:
            linelist = re.split(r" ", line)
            osspath = linelist[-1]
            
            osslist = osspath.strip("/").split("/")
            flowcell = osslist[-3]
            libid = osslist[-2]

            cmd = textwrap.dedent(
            """
            ossutil cp -r -u {osspath} {outdir}/{flowcell}/{libid}
            """.format(**locals()))

            print(flowcell, libid)
            print(cmd)


if __name__ == "__main__":
    
    args = sys.argv
    ossfile = args[1]
    outfile = args[2]
    outdir = args[3]

    cp(ossfile, outfile, outdir)

