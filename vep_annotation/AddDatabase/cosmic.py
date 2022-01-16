#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   cosmic.py
@Time    :   2021/11/13 11:40:43
@Author  :   lvmt 
@Version :   1.0
'''


from collections import defaultdict
from utils import utils



class COSMIC:
    
    def __init__(self, args):
        self.infile = args['infile']
        self.result = args['result']
        self.cosmic = args['cosmic']


    def get_cosmic_info(self):
        '''
        info = {
            'gene_chr_pos_ref_alt': line,
            'gene_chgvs': line
        }
        '''
        cosmic_info = defaultdict(dict)
        with utils.safe_open(self.cosmic, 'r') as fr:
            for line in fr:
                linelist = line.strip('\n').split('\t')
                gene = linelist[0]
                chgvs = linelist[2]
                _chr = linelist[4]
                pos = linelist[5]
                ref = linelist[6]
                alt = linelist[7]

                #cosmic为vcf格式的pos,ref,alt，更新为maf格式
                pos, ref, alt = utils.get_pos_ref_alt_from_vcf(pos, ref, alt)

                key1 = '{gene}_{chgvs}'.format(**locals())
                key2 = '{gene}_{_chr}_{pos}_{ref}_{alt}'.format(**locals())

                cosmic_info[key1] = line.strip('\n')
                cosmic_info[key2] = line.strip('\n')

        return cosmic_info

    
    def start(self):
        # print('>>>cosmic分析中')
        cosmic_info = self.get_cosmic_info()

        with open(self.infile, 'r') as fr, open(self.result, 'w') as fw:
            for line in fr:
                if line.startswith('Scale'):
                    head = line
                    fw.write(line)
                    head_index = utils.get_head_index(head)
                    continue

                linelist = line.strip('\n').split('\t')
                gene = linelist[head_index['hugo_symbol']]
                _chr = linelist[head_index['chromosome']]
                pos = linelist[head_index['start_position']]  #清楚为什么使用start!!
                ref = linelist[head_index['reference_allele']]
                alt = linelist[head_index['allele']]
                chgvs = linelist[head_index['hgvsc']]

                key1 = '{gene}_{_chr}_{pos}_{ref}_{alt}'.format(**locals())
                key2 = '{gene}_{chgvs}'.format(**locals())

                cosmic_tmp = ''

                tmp1 = cosmic_info.get(key1,'')
                tmp2 = cosmic_info.get(key2, '')

                if tmp1 == tmp2 and tmp1:
                    cosmic_linelist = tmp1.split('\t')
                    if cosmic_linelist[11] == '-':
                        cosmic_tmp = 'the mutation {0} has been exclude from the website.:{1}'.\
                                    format(cosmic_linelist[1], cosmic_linelist[10])
                    else:
                        cosmic_tmp = '{0}:{1};{2}'.format(cosmic_linelist[1], 
                                                    cosmic_linelist[10], 
                                                    cosmic_linelist[11])

                elif tmp1 and not tmp2:
                    cosmic_linelist = tmp1.split('\t')
                    if cosmic_linelist[11] == '-':
                        cosmic_tmp = 'the mutation {0} has beed excluded from \
                                      the website and cosmic gene or chgvs diff {1}_{2}:{3}'.format(
                                                                            cosmic_linelist[1],
                                                                            cosmic_linelist[0],
                                                                            cosmic_linelist[2],
                                                                            cosmic_linelist[10])
                                                                                                    
                    else:
                        cosmic_tmp = 'Cosmic gene or cHGVS diff {0}_{1}: {2}:{3};{4}'.format(
                                                                            cosmic_linelist[0], # gene
                                                                            cosmic_linelist[2], # hgvs
                                                                            cosmic_linelist[1], # cosm
                                                                            cosmic_linelist[10], # 1
                                                                            cosmic_linelist[11]) # large
                        
                elif tmp2 and not tmp1:
                    cosmic_linelist = tmp2.split('\t')
                    if cosmic_linelist[11] == '-':
                        cosmic_tmp = 'The mutation {0} has beed exclude from the website and \
                                      cosmic pos or alt diff {1}:{2} {3}/{4}: {5}'.format(
                                                                        cosmic_linelist[1],
                                                                        cosmic_linelist[4], # chr
                                                                        cosmic_linelist[5],
                                                                        cosmic_linelist[6],
                                                                        cosmic_linelist[7],
                                                                        cosmic_linelist[10])
                    else:
                        cosmic_tmp = 'Cosmic Pos or alt diff {0}:{1} {2}/{3}:{4}:{5};{6}'.format(
                                                                        cosmic_linelist[4], # chr
                                                                        cosmic_linelist[5],
                                                                        cosmic_linelist[6],
                                                                        cosmic_linelist[7],
                                                                        cosmic_linelist[1],
                                                                        cosmic_linelist[10],
                                                                        cosmic_linelist[11])
                
                else:
                    cosmic_tmp = '*'

                # 更细cosmic字段
                linelist[head_index['cosmic']] = cosmic_tmp
                fw.write('{0}\n'.format('\t'.join(linelist)))





