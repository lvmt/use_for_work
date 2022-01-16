#!/usr/bin/env python
#-*- coding:utf-8 -*-
'''
@Author: lvmengting
@Date: 2021-11-14 14:03:52
@Last Modified by:   lvmengting  
@Last Modified time: 2021-11-14 14:03:52
'''


from utils import utils

from logicalcheck.end_exon_check import EndExonCheck
from logicalcheck.gene_extend_check import GeneExtendCheck
from logicalcheck.splice_affect_check import SpliceAffectCheck
from logicalcheck.get_newfunction import NewFunction



class LogicalCheck:
    '''
    某些人工逻辑判断
    '''
    
    def parser_check(self, parser):
        parser.add_argument(
            '--infile', help='输入文件'
        )
        parser.add_argument(
            '--result', help='结果文件'
        )
        parser.add_argument(
            '--end_exon_config', help='末端exon配置文件'
        )
        parser.add_argument(
            '--gene_extend_confg', help='基因延长配置文件'
        )
        parser.add_argument(
            '--splice_affect_config', help='影响剪切配置文件'
        )
        
        parser.set_defaults(func=self.start)
        
        
    def start(self, **args):
        with open(args['infile'], 'r') as fr, open(args['result'], 'w') as fw:
            for line in fr:
                if line.startswith('Scale'):
                    fw.write(line)
                    head = line
                    head_index = utils.get_head_index(head)
                    continue

                linelist = line.strip('\n').split('\t')
                func = linelist[head_index['bgi_function']]
                tran = linelist[head_index['transcript_id']]
                start = linelist[head_index['start_position']]
                end = linelist[head_index['end_position']]
                phgvs = linelist[head_index['hgvsp_short']]
                chgvs = linelist[head_index['hgvsc']]
                flank = linelist[head_index['flank']]
                strand = linelist[head_index['strand']]
                gene = linelist[head_index['hugo_symbol']]
                exon_id = linelist[head_index['exon']]


                end_exon_tag = EndExonCheck(args, tran, start, end, func).start()
                gene_extend_tag = GeneExtendCheck(args, tran, func, phgvs).start()
                splice_affect_tag = SpliceAffectCheck(args, func, chgvs, flank, strand).start()
                newfun_tag = NewFunction(func, gene, exon_id, end_exon_tag, gene_extend_tag, splice_affect_tag).start()


                linelist[head_index['bgi_end_exon_check']] = end_exon_tag
                linelist[head_index['bgi_gene_extend_check']] = gene_extend_tag
                linelist[head_index['bgi_splice_affect_check']] = splice_affect_tag
                linelist[head_index['bgi_newfunction']] = newfun_tag
                # print(linelist)
                fw.write('{}\n'.format('\t'.join(linelist)))
 