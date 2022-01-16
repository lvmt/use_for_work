#!/usr/bin/env python
#-*- coding:utf-8 -*-
'''
@Author: lvmengting
@Date: 2021-11-14 13:53:40
@Last Modified by:   lvmengting  
@Last Modified time: 2021-11-14 13:53:40
'''

from utils import utils



class FilterReads:
    '''变异支持数过滤
    支持数>=3，并且正负链各一条
    '''
    def __init__(self, args):
        self.infile = args['infile']
        self.suffix = args['result_suffix']
        self.reads = args['reads']

    
    def pass_read_threshold(self, case_read, case_read_pos, case_read_neg):
        '''通过阈值条件
        '''
        if int(case_read) >= self.reads:
            if case_read_pos in ('*', '-'):
                return True
            elif int(case_read_pos) >= 1 and int(case_read_neg) >= 1:
                return True
            
        return False
    
    
    def start(self):
        '''
        '''
        with open(self.infile, 'r') as fr, open(self.suffix + '_pass', 'w') as fw_pass,\
            open(self.suffix + '_fail', 'w') as fw_fail:

            for line in fr:
                if line.startswith('Scale'):
                    head = line 
                    fw_pass.write(line)
                    fw_fail.write('{}\tfilter_reason\n'.format(line.strip('\n')))
                    head_index = utils.get_head_index(head)
                    continue

                linelist = line.strip('\n').split('\t')
                case_read = linelist[head_index['case_var_readsnum']]
                case_read_pos = linelist[head_index['case_var_positive_readsnum']]
                case_read_neg = linelist[head_index['case_var_negative_readsnum']]

                if self.pass_read_threshold(case_read, case_read_pos, case_read_neg):
                    fw_pass.write(line)
                else:
                    fw_fail.write('{}\tReadNum\n'.format(line.strip('\n')))

'''
fw_pass: 对应BGI的 1_snv_idl.anno.tsv
fw_fail: 对应BGI的 n_snv_idl.anno.tsv

read.list：来自输入文件= chr start
list: 来自于fw_pass, chr start end muttype gtype genotype # BGI和vep的坐标体系不一样,不知道对打分是否存在影响
'''
    
    