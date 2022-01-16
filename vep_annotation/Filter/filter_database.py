#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   filter_database.py
@Time    :   2022/01/05 11:04:05
@Author  :   lvmt 
@Version :   1.0
'''

## 要求：全部的库文件必须是vcf格式


# here put the import lib
from utils import utils
import os  


class FilterDataBase:

    def __init__(self, args):
        self.database_config = args['database_config']
        self.infile = args['infile']
        self.suffix = args['result_suffix']

    
    def get_database_info(self):
        '''将全部的database汇总成一个巨大哈希
        注意，此处的配置文件是vcf格式, vep为MAF格式
        在构建chr_pos_ref_alt，注意二者差异
        '''
        database_info = {}
        with open(self.database_config, 'r') as fr:
            for line in fr:
                linelist = line.strip('\n').split('\t')
                database_name = linelist[0]
                database_file = os.path.join(os.path.dirname(self.database_config), linelist[1])
                print(database_file)

                with open(database_file, 'r') as fr_database:
                    for line in fr_database:
                        if line.startswith('#'):
                            continue
                        linelist = line.strip('\n').split('\t')  #
                        _chr = linelist[0]
                        pos = linelist[1]
                        ref = linelist[2]
                        alt = linelist[3]
                        # 修改坐标体系为MAF格式，匹配vep注释结果
                        pos, ref, alt = utils.get_pos_ref_alt_from_vcf(pos, ref, alt)

                        key = '{_chr}_{pos}_{ref}_{alt}'.format(**locals())

                        database_info[key] = database_name 
        
        return database_info


    def start(self):
        database_info = self.get_database_info()

        with open(self.infile, 'r') as fr, open(self.suffix + '_pass', 'w') as fw_pass, \
            open(self.suffix + '_fail', 'w') as fw_fail:
            for line in fr:
                if line.startswith('Scale'):
                    fw_pass.write(line)
                    fw_fail.write('{}\tfail_reason\n'.format(line.strip('\n')))
                    head = line
                    head_index = utils.get_head_index(head)
                    continue
                linelist = line.strip('\n').split('\t')
                _chr = linelist[head_index['chromosome']]
                pos = linelist[head_index['start_position']]
                ref = linelist[head_index['reference_allele']]
                alt = linelist[head_index['allele']]

                key = '{_chr}_{pos}_{ref}_{alt}'.format(**locals())
                if key in database_info:
                    fw_fail.write('{}\t{}\n'.format(line.strip('\n'), database_info[key]))
                else:
                    fw_pass.write(line)

       
