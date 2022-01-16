#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   uniport.py
@Time    :   2021/11/19 11:40:43
@Author  :   lvmt 
@Version :   1.0
'''


import re
from collections import defaultdict
from utils import utils



class UNIPORT:
    '''更新输入文件中的三个关于uniport字段
    '''
    
    def __init__(self, args):
        self.infile = args['infile']
        self.resullt = args['result']
        self.uniport = args['uniport']


    def get_uniport_info(self):
        '''
        **uniport需要进行格式处理，包含5列信息
        获取uniport系列信息
        info = {
            'ABCB1'_51: {
                'length': 1280[51-357],
                'feature_key': DOMAIN,
                description: ABC transmembrane type-1 1.
            },

            TP53_283: {

            }
        }
        '''
        uniport_info = defaultdict(dict)
        uniport_gene_length = defaultdict()  # 记录uniport基因的长度

        with utils.safe_open(self.uniport, 'r') as fr:
            for line in fr:
                linelist = line.strip('\n').split('\t')
                gene = linelist[0]
                pos = linelist[1]
                length = linelist[2]
                feature = linelist[3]
                descript = linelist[4]

                key = '{gene}_{pos}'.format(**locals())
                uniport_info[key] = {
                    'length': length,
                    'feature_key': feature,
                    'description': descript
                }

                uniport_gene_length[gene] = length.split('[')[0]

        return uniport_info, uniport_gene_length
    
    
    def start(self):
        '''全部的最终接口均为start函数
        '''
        uniport_info, uniport_gene_length = self.get_uniport_info()

        with open(self.infile, 'r') as fr, open(self.resullt, 'w') as fw:
            for line in fr:
                if line.startswith('Scale'):
                    fw.write(line)
                    head = line
                    head_index = utils.get_head_index(head)
                    continue

                linelist = line.strip('\n').split('\t')
                gene = linelist[head_index['hugo_symbol']]
                phgvs = linelist[head_index['hgvsp_short']]

                #获取氨基酸发生突变的位置
                if re.search(r'(\d+)', phgvs):
                    phgvs_pos = re.search(r'(\d+)', phgvs).group(1)
                else:
                    phgvs_pos = '*'

                key = '{gene}_{phgvs_pos}'.format(**locals())

                linelist[head_index['bgi_uniport_position(s)']] = uniport_gene_length.get(gene, '*')

                if uniport_info.get(key):
                    linelist[head_index['bgi_uniport_position(s)']] = uniport_info[key]['length'] or '*'
                    linelist[head_index['bgi_uniport_feature_key']] = uniport_info[key]['feature_key'] or '*'
                    linelist[head_index['bgi_uniport_description']] = uniport_info[key]['description'] or '*'
                


                fw.write('{0}\n'.format('\t'.join(linelist)))

