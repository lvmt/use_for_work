#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   func_region.py
@Time    :   2021/11/29 16:25:50
@Author  :   lvmt 
@Version :   1.0
'''

# here put the import lib
from utils import utils



class FuncRegion:
    '''变异位于对应转录本的哪个CDS区域
    EX1 - CDS1
    or
    EX2 - CDS1
    '''
    def __init__(self, args):
        self.infile = args['infile']
        self.result = args['result']
        self.func_region_relation = args['func_region_relation']


    def get_func_relation_info(self):
        '''获取对应关系
        {
            'gene_trans_ex': cds
        }
        ''' 
        func_relation_info = {}
        with open(self.func_region_relation, 'r') as fr:
            for line in fr:
                linelist = line.strip('\n').split('\t')
                gene = linelist[0]
                tran = linelist[1]
                exon = linelist[2]
                cds = linelist[3]
                func_relation_info['{gene}_{tran}_{exon}'.format(**locals())] = cds 
        
        return func_relation_info


    def start(self):
        func_relation_info = self.get_func_relation_info()

        with open(self.infile, 'r') as fr, open(self.result, 'w') as fw:
            for line in fr:
                if line.startswith('Scale'):
                    head = line
                    fw.write(line)
                    head_index = utils.get_head_index(head)
                    continue
                linelist = line.strip('\n').split('\t')
                gene = linelist[head_index['hugo_symbol']]
                tran = linelist[head_index['transcript_id']]
                exon = linelist[head_index['exon']]

                key = '{gene}_{tran}_{exon}'.format(**locals())
                linelist[head_index['funcregion']] = func_relation_info.get(key, 'wrong')

                fw.write('{}\n'.format('\t'.join(linelist)))

