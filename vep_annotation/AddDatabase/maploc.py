#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   maploc.py
@Time    :   2021/11/29 15:33:54
@Author  :   lvmt 
@Version :   1.0
'''

# here put the import lib
from os import stat
from posix import listdir
import subprocess
from collections import defaultdict
from typing import ContextManager
from utils import utils


class MapLoc:
    '''变异位于哪个区间内
    '''

    def __init__(self, args):
        self.infile = args['infile']
        self.result = args['result']
        self.maploc = args['maploc']


    def get_maploc_info(self):
        '''{
            chr: [(start, end, loc), ()]
        }
        start+1，坐标体系是0base
        '''
        maploc_info = defaultdict(list)

        with open(self.maploc, 'r') as fr:
            for line in fr:
                linelist = line.strip('\n').split('\t')
                _chr = linelist[0]
                start = int(linelist[1])
                end = int(linelist[2])
                location = linelist[3]

                maploc_info[_chr].append((start + 1, end, location))

        return maploc_info

    
    def is_overlap(self, left1, right1, left2, right2):
        '''验证region1和region2是否存在交集
        '''
        if max(left1, left2) <= min(right1, right2):  # snv
            return True
        return False


    def get_variant_maploc(self, _chr, start, end, maploc_info):
        '''
        chr_start_end 
        目前先按照start,end进行注释，暂未发现跨2个区域的
        后续这个点待优化
        '''
        location_list = maploc_info[_chr]
        for (left2, right2, location) in location_list: 
            if self.is_overlap(start, end, left2, right2):
                return location

    
    def start(self):
        maploc_info = self.get_maploc_info()
        with open(self.infile, 'r') as fr, open(self.result, 'w') as fw:
            for line in fr:
                if line.startswith('Scale'):
                    fw.write(line)
                    head = line
                    head_index = utils.get_head_index(head)
                    continue
                linelist = line.strip('\n').split('\t')
                _chr = linelist[head_index['chromosome']]
                start = int(linelist[head_index['start_position']])
                end = int(linelist[head_index['end_position']])

                map_location = self.get_variant_maploc(_chr, start, end, maploc_info)
                linelist[head_index['maploc']] = map_location
                # print(linelist)
                fw.write('{}\n'.format('\t'.join(linelist)))








