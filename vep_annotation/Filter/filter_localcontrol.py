#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   filter_localcontrol.py
@Time    :   2021/12/03 13:59:33
@Author  :   lvmt 
@Version :   1.0
'''

# 要求：snv indel 库文件必须是标准vcf格式


# here put the import lib
from typing import List
from utils import utils
import re


    
class FilterLocalControl:
    '''
    基于本地位点频率数据库进行过滤
    保留：case_var_freq > 本地数据库频率
    '''

    def __init__(self, args):
        self.infile = args['infile']
        self.suffix = args['result_suffix']
        self.local_snv = args['local_snv']
        self.loval_indel = args['local_indel']


    def vcf2dict(self, vcf):
        info = {}
        with utils.safe_open(vcf, 'r') as fr:
            for line in fr:
                if line.startswith('#'):
                    continue
                linelist = line.strip('\n').split('\t')
                _chr = linelist[0]
                pos = linelist[1]
                ref = linelist[3]
                alt = linelist[4]
                tag = linelist[-1] # pop=0.284101599247413;fre=5.15384292456624

                freq = float(re.search(r'fre=(.*)', tag).group(1))
                pos, ref, alt = utils.get_pos_ref_alt_from_vcf(pos, ref, alt)
                info['{_chr}_{pos}_{ref}_{alt}'.format(**locals())] = [freq, tag]

        return info


    def get_localcontrol_info(self):
        localcontrol_info = {}
        snv_info = self.vcf2dict(self.local_snv)
        indel_info = self.vcf2dict(self.loval_indel)
        localcontrol_info.update(snv_info)
        localcontrol_info.update(indel_info)

        return localcontrol_info


    def start(self):
        localcontrol_info = self.get_localcontrol_info()

        with open(self.infile, 'r') as fr, open(self.suffix + '_pass', 'w') as fw_pass,\
            open(self.suffix + '_fail', 'w') as fw_fail:
            for line in fr:
                if line.startswith('Scale'):
                    head = line
                    fw_pass.write(line)
                    fw_fail.write(line)
                    head_index = utils.get_head_index(head)
                    continue
                
                linelist = line.strip('\n').split('\t')
                case_freq = float(linelist[head_index['case_var_freq']])
                _chr = linelist[head_index['chromosome']]
                pos = linelist[head_index['start_position']]  # 理解为什么选取start作为pos
                ref = linelist[head_index['reference_allele']]
                alt = linelist[head_index['allele']]

                key = '{_chr}_{pos}_{ref}_{alt}'.format(**locals())

                if key in localcontrol_info and case_freq > localcontrol_info[key][0]:
                    fw_pass.write('NOControl{0};{1}'.format(localcontrol_info[key][1], line))
                elif key in localcontrol_info and case_freq <= localcontrol_info[key][0]:
                    fw_fail.write('NOControl{0};{1}'.format(localcontrol_info[key][1], line))
                else:
                    fw_pass.write(line)

