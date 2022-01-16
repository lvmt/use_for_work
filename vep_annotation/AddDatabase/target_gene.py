#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   target_gene.py
@Time    :   2021/11/19 11:27:01
@Author  :   lvmt 
@Version :   1.0
'''

# here put the import lib
from collections import defaultdict
from utils import utils



class TargetGene:

    def __init__(self, args):
        self.infile = args['infile']
        self.result = args['result']
        self.gene_pos = args['target_gene_pos']
        self.gene_func = args['target_gene_func']
        self.gene_exon = args['target_gene_exon']
        

    def get_gene_pos_info(self):
        '''
        info = {
            'gene_p.E37K'：{
                'Tier II': '非小细胞肺癌;结直肠癌',
                'Tier I': '甲状腺癌'
            },

        }
        '''
        chgvs_info = defaultdict(dict)
        phgvs_info = defaultdict(dict)

        with open(self.gene_pos, 'r') as fr:
            for line in fr:
                if line.startswith("#"):
                    continue
                linelist = line.strip('').split('\t')
                gene = linelist[0]
                chgvs = linelist[1]
                phgvs = linelist[2]
                variant_class = linelist[3]
                caner_class = linelist[4]

                if chgvs:
                    key1 = '{gene}_{chgvs}'.format(**locals())
                    if not key1 in chgvs_info.keys():
                        chgvs_info[key1] = defaultdict(str)
                    chgvs_info[key1][variant_class] += caner_class + ';'

                if phgvs:
                    key2 = '{gene}_{phgvs}'.format(**locals())
                    if not key2 in phgvs_info.keys():
                        phgvs_info[key2] = defaultdict(str)
                    phgvs_info[key2][variant_class] += caner_class + ';'

        return chgvs_info, phgvs_info


    def get_gene_func_info(self):
        '''
        info = {
            'PTEN_frameshift': {
                'Tier II': '非小细胞肺癌;软组织肉瘤',
                'Tier I': '乳腺癌'
            },
            
        }
        '''
        func_info = defaultdict(dict)
        with open(self.gene_func, 'r') as fr:
            fr.readline()
            for line in fr:
                if line.startswith('#'):
                    continue
                linelist = line.strip('\n').split('\t')
                gene = linelist[0]
                muttype = linelist[3] # frameshift
                variant_class = linelist[1]
                cancer_class = linelist[2]
                key = '{gene}_{muttype}'.format(**locals())

                if not key in func_info.keys():
                    func_info[key] = defaultdict(str)
                    
                func_info[key][variant_class] += cancer_class + ';'

        return func_info


    def get_gene_exon_info(self):
        exon_info = defaultdict()
        with open(self.gene_exon, 'r') as fr:
            for line in fr:
                linelist = line.strip('').split('\t')
                gene = linelist[0]
                exon_info[gene] = 1
        
        return exon_info


    def start(self):
        func_info = self.get_gene_func_info()
        chgvs_info, phgvs_info = self.get_gene_pos_info()
        exon_info = self.get_gene_exon_info()

        with open(self.infile, 'r') as fr, open(self.result, 'w') as fw:
            for line in fr:
                if line.startswith('Scale'):
                    head = line
                    fw.write(line)
                    head_index = utils.get_head_index(line)
                    continue

                linelist = line.strip('\n').split('\t')
                gene = linelist[head_index['hugo_symbol']]
                chgvs = linelist[head_index['hgvsc']]
                phgvs = linelist[head_index['hgvsp_short']]
                func = linelist[head_index['bgi_function']]
                exon = linelist[head_index['exon']]

                key_func = '{gene}_{func}'.format(**locals())
                key_chgvs = '{gene}_{chgvs}'.format(**locals())
                key_phgvs = '{gene}_{phgvs}'.format(**locals())

                # print(key_func)
                #修改Target_gene字段信息
                tmp = ''
                if func_info.get(key_func):
                    for key,value in func_info[key_func].items():
                        tmp += '{key}:{value};'.format(**locals())
                elif phgvs_info.get(key_phgvs):
                    for key,value in phgvs_info[key_phgvs].items():
                        tmp += '{key}:{value};'.format(**locals())
                elif chgvs_info.get(key_chgvs):
                    for key,value in chgvs_info[key_chgvs].items():
                        tmp += '{key}:{value};'.format(**locals())
                elif exon_info.get(gene) and exon.startswith('EX'):
                    tmp = exon_info[gene]

                # tmp = tmp.replace(';;', ';')

                if tmp:
                    linelist[head_index['target_gene']] = 'YES({tmp})'.format(**locals()).replace(';;', ';')
                else:
                    linelist[head_index['target_gene']] = 'NO'

                fw.write('{}\n'.format('\t'.join(map(str, linelist))))



                

                
                





                

