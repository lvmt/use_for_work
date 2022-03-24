#!/usr/bin/env python3
#-*- coding:utf-8 -*-


'''
基于gene和转录本信息，从旧版本的gff文件中，提取指定转录本的gff文件
1）基于gene名字获取gene列信息

'''

import subprocess
import sys 
import re  


gff = sys.argv[1]
gene_file = sys.argv[2]
out_gff = sys.argv[3]

with open(gene_file, 'r') as fr, open(out_gff, 'w') as fw:
    for line in fr:
        linelist = line.strip('\n').split('\t')
        gene = linelist[0]
        tran = linelist[1]

        cmd_gene = subprocess.getoutput("less {gff}| grep 'Name={gene};' | awk 'if($3==\"gene\") print $0'".format(**locals()))
        cmd_tran = subprocess.getoutput("less {gff}| grep {tran} | awk '{if($3==\"mRNA\" || $3 == \"transcript\") print $0}".format(**locals()))
        rna_id = re.search(r'ID=(.*?);', cmd_tran).group(1)
        exon_cds = subprocess.getoutput('less {gtf} | grep {rna_id}').split('\n')

        fw.write(cmd_gene + '\n')
        fw.write(cmd_tran + '\n')
        for exon in exon_cds:
            fw.write(exon + '\n')

            
        



