#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   get_newfunction.py
@Time    :   2021/12/13 15:47:06
@Author  :   lvmt 
@Version :   1.0
'''

# here put the import lib
'''
解读Newfunction字段的确定
是根据bgifunction加上以下3个字段的结果进行更新的
 - EndExonCheck
 - GeneExtendChek
 - SpliceAffectCheck

'''



from threading import local


class NewFunction:

    # 按照行的格式进行处理
    
    def __init__(self, func, gene, exon_id, end_exon, gene_extend, splice_affect):
        self.func = func
        self.gene = gene 
        self.exon_id = exon_id
        self.end_exon = end_exon
        self.gene_extend = gene_extend
        self.splice_affect = splice_affect


    def handle_frameshift(self):
        # gene_extend has the best priority
        if self.gene_extend == 'Y':
            return 'frameshift-L'
        elif self.end_exon == 'Y':
            return 'frameshift-N'
        else:
            return 'frameshift'


    def handle_nonsense(self):
        if self.end_exon == 'Y':
            return 'nonsense-N'
        return 'nonsense'


    def handle_splice(self):
        if self.splice_affect == 'Y' and self.end_exon == 'N':
            return '{self.func}-Y'.format(**locals())
        elif self.splice_affect == 'Y' and self.end_exon == 'Y':
            return '{self.func}-N'.format(**locals())
        return self.func


    def handle_span(self):
        if self.splice_affect == 'Y' and self.end_exon == 'N':
            return '{self.func}-Y'.format(**locals())
        elif self.splice_affect == 'Y' and self.end_exon == 'Y':
            return '{self.func}-N'.format(**locals())
        return self.func


    def get_newfunction(self):
        func_dict = {
            'frameshift': self.handle_frameshift(),
            'nonsense': self.handle_nonsense(),
            'splice-3': self.handle_splice(),
            'splice-5': self.handle_splice(),
            'span': self.handle_span()
        }
    
        return func_dict.get(self.func, self.func)

    
    def handle_MET(self, newfunc):
        '''如果是MET基因,对于起特定区域的newfunc需要进行重新矫正
        '''
        target_met_set = set(('IVS13', 'IVS13-EX14', 'EX14-IVS13', 'EX14', 'EX14-IVS14', 'IVS14-EX14', 'IVS14'))
        exon_id_set = set(self.exon_id.split('-'))

        if 'Y' in newfunc and (target_met_set & exon_id_set):
            newfunc = newfunc.replace('Y', 'M')
        
        return newfunc


    def start(self):
        newfunc = self.get_newfunction()
        
        if self.gene == 'MET':
            return self.handle_MET(newfunc)

        return newfunc
            
