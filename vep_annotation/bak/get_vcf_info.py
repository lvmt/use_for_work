#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   get_vcf_info.py
@Time    :   2021/11/12 14:53:43
@Author  :   lvmt 
@Version :   1.0
'''

# here put the import lib
# 处理vcf文件，得到变异的支持数信息
# snv和indel获取变异支持数的时候，读取的文件是不一样的
# 需要适配单双组织



import re
from collections import defaultdict

import utils



class HandleVcf:



    def __init__(self, vcf, filetype):
        self.vcf = vcf
        self.filetype = filetype  # snv or indel
        self.case_pattern = re.compile(r'DP4=(\d+),(\d+),(\d+),(\d+);?')
        self.control_pattern = re.compile(r'CtrlDP4=(\d+),(\d+),(\d+),(\d+);?')

    
    def get_snv_info(self):
        '''处理snv vcf文件
        读取INFO字段,获取需要的信息
        vcf_info = {
            chr_pos_ref/alt: [(case), (control)]
        }
        key的形式和vep注释结果中的upload_variation字段一致
        '''
        vcf_info = defaultdict(dict)
        with utils.safe_open(self.vcf, 'r') as fr:
            for line in fr:
                if line.startswith('##'):
                    continue
                elif line.startswith('#'):
                    head = line.strip('')
                    head_index = utils.get_head_index(head)
                    continue

                linelist = line.strip('').split('\t')
                _chr = linelist[head_index['#chrom']]
                pos = linelist[head_index['pos']]
                ref = linelist[head_index['ref']]
                alt = linelist[head_index['alt']]
                key = '{_chr}_{pos}_{ref}/{alt}'.format(**locals())

                info = linelist[head_index['info']]
                case = re.search(self.case_pattern, info)
                control = re.search(self.control_pattern, info)
                if case:
                    case_ref_positive = int(case.group(1))
                    case_ref_negative = int(case.group(2))
                    case_alt_positive = int(case.group(3))  # read2 is alt read
                    case_alt_negative = int(case.group(4))
                    case_ref = case_ref_positive + case_ref_negative
                    case_alt = case_alt_positive + case_alt_negative
                    case_pos_depth = case_ref + case_alt
                    case_var_freq = (case_alt / case_pos_depth) * 100
                    vcf_info[key].update({
                        'Case_ref_readsNum': case_ref,
                        'Case_var_readsNum': case_alt,
                        'Case_var_Positive_readsNum': case_alt_positive,
                        'Case_var_Negative_readsNum': case_alt_negative,
                        'Case_pos_dep': case_pos_depth,
                        'Case_var_freq': case_var_freq
                    })
                if control:
                    control_ref_positve = int(control.group(1))
                    control_ref_negative = int(control.group(2))
                    control_alt_positve = int(control.group(3))
                    control_alt_negative = int(control.group(4))
                    control_ref = control_ref_positve + control_ref_negative
                    control_alt = control_alt_positve + control_alt_negative
                    control_pos_depth = control_ref + control_alt
                    vcf_info[key].update({
                        'Ctrl_ref_readsNum': control_ref,
                        'Ctrl_var_readsNum': control_alt,
                        'Ctrl_var_Positive_readsNum': control_alt_positve,
                        'Ctrl_var_Negative_readsNum': control_alt_negative,
                        'Ctrl_pos_dep': control_pos_depth
                    })

        return vcf_info


    def get_indel_info(self):
        '''处理indel vcf文件
        处理indel的vcf文件时，为了和注释的pos,以及ref和alt对应，需要对vcf做处理
        {
            'chr_pos_ref/alt': [(case), (control)]
        }
        '''
        with utils.safe_open(self.vcf, 'r') as fr:
            for line in fr:
                if line.startswith('##'):
                    pass
                elif line.startswith('#'):
                    head = line.strip('')
                    head_index = utils.get_head_index(head)
                    continue

                linelist = line.strip('').split('\t')
                _chr = linelist[head_index['#chrom']]
                pos = linelist[head_index['pos']]
                ref = linelist[head_index['ref']]
                alt = linelist[head_index['alt']]

                pos, ref, alt = utils.modify_pos_ref_alt(pos, ref, alt)
                if 'cancer' in head_index:
                    case = linelist[head_index['cancer']]
                if 'normal' in head_index:
                    control = linelist[head_index['normal']]


    def start(self):
        func_info = {
            'snv': self.get_snv_info,
            'indel': self.get_indel_info
        }

        freq_info = func_info[self.filetype]()

        return freq_info






