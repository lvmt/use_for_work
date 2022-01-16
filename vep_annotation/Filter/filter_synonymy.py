#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   filter_synonymy.py
@Time    :   2022/01/05 14:52:26
@Author  :   lvmt 
@Version :   1.0
'''

# here put the import lib
# 过滤掉同义突变--只过滤掉TMB标签为NoTMB的同义突变
from types import CoroutineType
from utils import utils


class FilterSynonymy:

    def __init__(self, args):
        self.infile = args['infile']
        self.suffix = args['result_suffix']

    
    def start(self):
        with open(self.infile, 'r') as fr, open(self.suffix + '_pass', 'w') as fw_pass,\
            open(self.suffix + '_fail', 'w') as fw_fail:
            for line in fr:
                if line.startswith('Scale'):
                    fw_pass.write(line)
                    fw_fail.write('{}\tfilter_reason\n'.format(line.strip('\n')))
                    head = line
                    head_index = utils.get_head_index(head)
                    continue
                linelist = line.strip('\n').split('\t')
                tmb_tag = linelist[head_index['tmb_type']]
                func = linelist[head_index['bgi_function']]

                if func == 'coding-synon' and tmb_tag == 'noTMB':
                    fw_fail.write('{}\tsynonymy\n'.format('\t'.join(linelist)))
                else:
                    fw_pass.write(line)
                    
                    