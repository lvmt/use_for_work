#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   end_exon_check.py
@Time    :   2021/12/13 15:41:45
@Author  :   lvmt 
@Version :   1.0
'''

# here put the import lib
'''
末端exon检查
check标签
 - frameshift
 - ninsense
 - splice-3
 - splice-5
 - span

其他标签输出Nan
判断是否延长的条件 ：
变异完全落在end_exon的end后面，如果是以overlap形式存在，不认为是end_exon

逻辑：
+ mut_start >= region_start
- ：mut_end <= region_end  
'''



from collections import defaultdict
from os import truncate
from posix import EX_NOHOST
import re


class EndExonCheck:
    
    def __init__(self, args, tran, start, end, func):
        self.end_exon_database = args['end_exon_config']
        self.tran = tran 
        self.mut_start = start
        self.mut_end = end 
        self.func = func


    def get_end_exon_info(self):
        end_exon_info = defaultdict()
        with open(self.end_exon_database, 'r') as fr:
            for line in fr:
                if line.startswith('#'):
                    continue
                linelist = line.strip('\n').split('\t')
                start = linelist[1]
                end = linelist[2]
                tran = linelist[4] 
                strand = linelist[5]
                if not tran in end_exon_info:
                    end_exon_info[tran] = [start, end, strand]

                if start < end_exon_info[tran][0]:
                    end_exon_info[tran][0] = start
                if end > end_exon_info[tran][1]:
                    end_exon_info[tran][1] = end 

        return end_exon_info


    def positive_end_exon(self, mut_start, mut_end, region_start, region_end):
        if mut_start >= region_start:
            return 'Y'
        return 'N'

    def negative_end_exon(self, mut_start, mut_end, region_start, region_end):
        if mut_end <= region_end:
            return 'Y'
        return 'N'


    def start(self):
        end_exon_info = self.get_end_exon_info()
        region_start, region_end, strand = end_exon_info[self.tran]

        target_func = ['frameshift', 'nonsense', 'splice-3', 'splice-5', 'span']
        
        if self.func in target_func and strand == '+':
            return self.positive_end_exon(self.mut_start, self.mut_end, region_start, region_end)
        elif self.func in target_func and strand == '-':
            return self.negative_end_exon(self.mut_start, self.mut_end, region_start, region_end)
        else:
            return 'Nan'
       


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='just a test')
    parser.add_argument('--end_exon_config', help='end_exon库文件')
    parser.add_argument('--infile', help='test file')

    args = vars(parser.parse_args())

    with open(args['infile'], 'r') as fr:
        for line in fr:
            if line.startswith('Scale'):
                continue
            linelist = line.strip('\n').split('\t')
            tran = linelist[37]
            start = linelist[29]
            end = linelist[30]
            func = linelist[9]

            print(tran, start, end, func, EndExonCheck(args, tran, start, end, func).start())
