#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   transver_function.py
@Time    :   2021/11/15 14:45:09
@Author  :   lvmt 
@Version :   1.0
'''

# here put the import lib
'''将vep原始注释功能转换为对应的BGICG功能
1. 先根据规则注释功能进行转换
2. 在根据功能yaml文件,进行替换
'''


'''
注释结果中2类特殊function需要优先处理
splice_region_variant: 遇到这个func,直接从func_list中del
protein_altering_variant：第一优先级,需要根据氨基酸的phgvs进行进一步判断
'''


import yaml
import re
# from simplify import simplify_vep_annotation




class TransverFunction:

    def __init__(self, **args):
        self.vep_function_yaml = args['vep_function_yaml']


    def simplify_function(self, vep_function, tert, gene):
        '''
        优先处理tert 启动子区域问题
        vep的注释结果可能存在多个,需要简化成为一个
        '''
        #因为splice_region_variant是个陪衬,先清除掉
        vep_function = vep_function.replace('splice_region_variant', '')
        vep_function = vep_function.strip(',') # 去除多余的分隔符 &

        func_list = vep_function.split(',')
        if len(func_list) == 1 and gene == 'TERT' and tert.startswith('NM'):
            return 'promter'
        elif len(func_list) == 1:
            return vep_function
        elif len(func_list) > 1:  # 及时func_list=3, 也不一定是span
            return self.get_the_better_func(func_list)

    
    def get_the_better_func(self, func_list):
        '''
        1.有时候2个注释也会是span，因为跨过了不同的编码区域
        2.根据function的优先级,返回优先级较高的func信息
        '''
        region_set = set()
        func_region_info = yaml.load(open(self.vep_function_yaml))['FuncRegion']
        func_priority = yaml.load(open(self.vep_function_yaml))['Priority']

        for func in func_list:
            region_set.add(func_region_info[func])

        if len(region_set) > 1:
            return 'span'
        else:
            if func_priority[func_list[0]] < func_priority[func_list[1]]:
                return func_list[0]
            else:
                return func_list[1]


    def vep2bgi(self, simplify_vep_func, hgvsc, hgvsp, ref, alt):
        '''将vep的func转化为bgi的func
        优先处理protein_altering_variant
        优先处理coding_sequence_variant
        对于span需要进行再次矫正
        '''
        vep2bgi_info = yaml.load(open(self.vep_function_yaml))['vep2BGI']

        if simplify_vep_func == 'protein_altering_variant':
            return self.handle_protein_altering_variant_from_hgvsp(hgvsp, ref, alt)
        elif simplify_vep_func == 'coding_sequence_variant':
            return 'coding_sequence_variant'
            # return self.handle_coding_sequence_variant_from_chgvs_phgvs(hgvsc, hgvsp)
        elif simplify_vep_func == 'span':
            return self.correct_span(hgvsp, ref, alt)
        else:
            return vep2bgi_info[simplify_vep_func]


    def handle_protein_altering_variant_from_hgvsp(self, hgvsp, ref, alt):
        '''根据phgvs判断bgi的function结果
        现根据phgvs判断:stop_gain,frameshift
        在根据ref，alt判断: 
        cds-ins, cds-del,cds-indel  
        '''
        if 'fs' in hgvsp:
            return 'frameshift'
        elif hgvsp.endswith('Ter'):
            # stop_gained
            return 'nonsense'
        elif 'delins' in hgvsp:  #还需进一步判断,解读逻辑白帅帅
            if len(ref) > len(alt):
                return 'cds-del'
            elif len(ref) < len(alt):
                return 'cds-ins'
            else:
                return 'null'
        elif 'del' in hgvsp:
            return 'cds-del'
        elif 'ins' in hgvsp:
            return 'cds-ins'
        else:
            return 'null'


    def handle_coding_sequence_variant_from_chgvs_phgvs(self, hgvsc, hgvsp, ref, alt):
        '''
        神马玩意
        竟然单独出现了coding_sequence_variant
        决策就是根据chgvs和phgvs进行判断
        若存在phgvs,则利用handle_proten***这个函数进行判断
        否则按照chgvs进行判断
        c.3727_3727+1delinsAC 
        好像还有更复杂的格式问题，需要解决：20211117
        这个也有部分问题待解决
        '''
        if not hgvsp == '-':
            return self.handle_protein_altering_variant_from_hgvsp(hgvsp, ref, alt)
        else:
            # 内含子区域或者剪切区
            flag, num = re.search(r'([+-])(\d+)[A-Z]', hgvsc).groups()
            if int(num) > 2:
                return 'intron'
            elif flag == '+':
                return 'splice-5'
            elif flag == '-':
                return 'splice-3'
            else:
                return 'null'

            
    def correct_span(self, hgvsp, ref, alt):
        '''
        vep中出现一类特殊案例
        splice_acceptor_variant&frameshift_variant&intron_variant
        基于上述的逻辑判断，最终的功能会注释为span，但是在回查注释结果时，phgvs是由结果的
        因此对于注释为span的类型，我们需要进行再次矫正
        '''
        if hgvsp == '-':
            return 'span'
        else:
            return self.handle_protein_altering_variant_from_hgvsp(hgvsp, ref, alt)