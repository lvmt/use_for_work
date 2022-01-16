#!/usr/bin/env python
#-*- coding:utf-8 -*-
'''
@Author: lvmengting
@Date: 2021-11-14 13:51:49
@Last Modified by:   lvmengting  
@Last Modified time: 2021-11-14 13:51:49
'''


from utils import utils



class FilterFreqDatabase:
    '''基于频率数据库注释
    '''
    def __init__(self, args):
        self.infile = args['infile']
        self.suffix = args['result_suffix']
        self.freq = args['freq']
        self.freq_database_ratio = args['freq_database_ratio']

    
    def get_freq_item(self):
        '''基于哪些freq数据库进行注释
        至少一个库
        '''
        freq_list = [
            'GMAF',
            'EAS_AF',
            'ESP6500',
            'gnomAD_AF',
            'gnomAD_AFR_AF'
        ]
        
        return freq_list


    def judge_freq_tag(self, freq_status_list):

        '''在freq_status_list,有多少个满足情况保留
        假设self.freq的阈值为1，大于1的status为fail，小于status为pass
        freq_status_list： [pass, pass, fail, pass, fail]
        threshold为列表中满足条件的数目,假设为3,则返回True
        '''
        # 
        threshold = len(freq_status_list) * self.freq_database_ratio

        if freq_status_list.count('pass') >= threshold:
            return True
        return False

    
    def start(self):
        freq_item = self.get_freq_item()

        with open(self.infile, 'r') as fr, open(self.suffix + '.filter_freq_pass', 'w') as fw_pass,\
             open(self.suffix + '.filter_freq_fail', 'w') as fw_fail:

            for line in fr:
                freq_status_list = []  # 
                if line.startswith('Scale'):
                    fw_pass.write(line)
                    fw_fail.write(line)
                    head = line
                    head_index = utils.get_head_index(head)
                    continue

                linelist = line.strip('').split('\t')
                for item in freq_item:  # item = gnomAD_SAS_AF
                    if not linelist[head_index[item.lower()]] == '-' and\
                         float(linelist[head_index[item.lower()]]) < self.freq:

                        freq_status_list.append('pass')
                    else:
                        freq_status_list.append('fail')

                if self.judge_freq_tag(freq_status_list):
                    fw_pass.write(line)
                else:
                    fw_fail.write(line)

                # 清空，列表为可更改对象
                freq_status_list.clear()
    