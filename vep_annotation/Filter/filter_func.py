#!/usr/bin/env python
#-*- coding:utf-8 -*-
'''
@Author: lvmengting
@Date: 2021-11-14 13:52:44
@Last Modified by:   lvmengting  
@Last Modified time: 2021-11-14 13:52:44
'''


from numpy import save
import yaml
import re
from utils import utils


class FilterFunc:
    '''功能过滤
    初步保留的功能列表见：function.yaml
    '''

    def __init__(self, **args):
        self.func_config = args['func_config']
        self.infile = args['infile']
        self.suffix = args['result_suffix']

    
    def get_save_function_list(self):
        '''
        需要保留的功能info
        '''
        func_list = yaml.load(open(self.func_config))
        return func_list
    
    
    def save_span(self, chgvs, bgi_func, vep_func):
        '''只保留特定类型的span
        只要其中一端位于保留区域即可以保留（基于chgvs）
        c.*100/c.-100 直接不要
        c.100+x (x<=2保留, 其他不要)
        20220324：某些特例NM_005243.3:c.1932-21_*9del,
        针对特例，需要利用vep_func，判断其中是否存在coding_sequence_variant
        '''
        if bgi_func == 'span' and chgvs == '-':  # 有些span可能没有chgvs注释,多位于UTR区域
            return False

        if bgi_func == 'span' and 'coding_sequence_variant' in vep_func:
            return True

        splice_p = re.compile(r'(\d+)[+-][12]($|[a-zA-Z])')  # 满足条件的+-1/2; 100-2或者100-2delins
        exon_p = re.compile(r'(\d+)($|[a-zA-Z])')  # 满足条件的exon; 100或者100delins

        if bgi_func == 'span':
            save_tag = []  # 记录span两端是否位于保留区域
            chgvs = chgvs.split('c.')[-1].split('_')

            for item in chgvs:
                if item.startswith('-') or item.startswith('*'):  # UTR区域,不参与判断
                    continue
                if ('-' in item or '+' in item) and re.findall(splice_p, item):
                    save_tag.append(True)
                elif (not '-' in item) and (not '+' in item) and re.findall(exon_p, item):
                    save_tag.append(True)
                else:
                    save_tag.append(False)

            if any(save_tag): # 任何一端位于保留区域,均保留下来
                return True
            else:
                return False

        return True
        
    
    def start(self):
        func_list = self.get_save_function_list()

        with utils.safe_open(self.infile, 'r') as fr, open(self.suffix + '_pass', 'w') as fw_pass, \
            open(self.suffix + '_fail', 'w') as fw_fail:
            
            for line in fr:
                if line.startswith('Scale'):
                    fw_pass.write(line)
                    fw_fail.write(line)
                    head = line 
                    head_index = utils.get_head_index(head)
                    continue
                linelist = line.strip('').split('\t')
                bgi_func = linelist[head_index['bgi_function']]  # nonsense
                vep_func = linelist[head_index['consequence']]  # 原始vep注释结果
                gene = linelist[head_index['hugo_symbol']] # TERT
                chgvs = linelist[head_index['hgvsc']] # c.-146C>T

                # 针对8个tert启动子，构造特殊key
                tert_promter = '{gene}:{chgvs}'.format(**locals()) 
                if bgi_func in func_list and self.save_span(chgvs, bgi_func, vep_func) or tert_promter in func_list:
                    fw_pass.write(line)
                else:
                    fw_fail.write(line)



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='test')
    parser.add_argument('--func_config')
    parser.add_argument('--infile')
    parser.add_argument('--result_suffix')

    args = vars(parser.parse_args())
    FilterFunc(**args).start()


 