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
from collections import defaultdict


class OneItemList:
    '''
    定义只有一个元素的列表结果
    '''
    def __init__(self):
        self.items = []


    def push(self, item):
        if len(self.items) > 1:
            self.items.clear()
        self.items.append(item)
        return self


    def final_exon(self):
        return self.items[0]



class FuncRegion:
    '''变异位于对应转录本的哪个CDS区域
    EX1 - CDS1
    or
    EX2 - CDS1
    '''
    def __init__(self, args):
        self.infile = args['infile']
        self.result = args['result']
        self.func_region_relation = args['func_region_relation']  # dbref格式


    def get_final_exon_info(self):
        '''
        final_exon_info = {'tran': ['exon100']}
        '''
        final_exon_info = {}

        with open(self.func_region_relation, 'r') as fr:
            for line in fr:
                linelist = line.strip('\n').split('\t')
                gene = linelist[3]
                tran = linelist[4]
                exon = linelist[6]
                if tran in final_exon_info:
                    final_exon_info[tran].push(exon)
                else:
                    final_exon_info[tran] = OneItemList()
                    final_exon_info[tran].push(exon)
        
        return final_exon_info


    def final_exon(self, final_exon_info, tran, exon):
        '''判断是否是最后一个exon'''
        if final_exon_info.get(tran):
            final_exon = final_exon_info[tran].final_exon()
            if exon == final_exon:
                return '{exon}E'.format(**locals())
        return exon


    def get_func_relation_info(self):
        '''获取对应关系
        {
            'trans_ex': cds
        }
        ''' 
        func_relation_info = {}
        with open(self.func_region_relation, 'r') as fr:
            for line in fr:
                linelist = line.strip('\n').split('\t')
                gene = linelist[3]
                tran = linelist[4]
                exon = linelist[6]
                cds = linelist[7]
                func_relation_info['{tran}_{exon}'.format(**locals())] = cds 
        
        return func_relation_info


    def start(self):
        final_exon_info = self.get_final_exon_info()
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

                key = '{tran}_{exon}'.format(**locals())
                # 更新Funcregion 字段
                linelist[head_index['funcregion']] = func_relation_info.get(key, 'Nan')
                # 更新最后一个exon的写法
                linelist[head_index['exon']] = self.final_exon(final_exon_info, tran, exon)

                fw.write('{}\n'.format('\t'.join(linelist)))




if __name__ == '__main__':
    # 测试OneItemList
    pass 


