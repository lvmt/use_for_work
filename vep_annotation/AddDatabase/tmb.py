#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   tmb.py
@Time    :   2021/11/29 16:27:49
@Author  :   lvmt 
@Version :   1.0
'''

# here put the import lib
'''
TMB计算的输入文件需要进行以下处理：
1. 功能区域过滤
2. 频率数据库过滤
3. read条数过滤

'''



import re 
import yaml
from utils import utils



class TMB:
    '''TMB_Type字段
    '''
    def __init__(self, args):
        self.infile = args['infile']
        self.result = args['result']  #增加TMB_Type字段 Driver, TMB, noTMB
        self.driver_list = args['driver_list']  # 驱动基因列表，库文件格式必须处理好
        self.driver_yaml = args['driver_yaml']  # 某些失传已久的delins基因
        self.tmb_freq = args['tmb_freq']  # 一般认为freq > 0.15
        self.chip_size = args['chip_size'] # 


    def get_driver_info(self):
        '''文件格式有点乱
        库文件中chgvs或者phgvs空值不能使用-指示
        第一列为gene，第二列为hgvs或者/，第三列为phgvs或者loss
        {
            'gene_chgvs': 1,
            'gene_phgvs': 1,
            'gene_loss': 1,
            'gene_/': 1
        }
        '''
        driver_info = {}
        with open(self.driver_list, 'r') as fr:
            for line in fr:
                linelist = line.strip('\n').split('\t')
                gene = linelist[0]
                chgvs = linelist[1]
                phgvs = linelist[2]

                chgvs_list = chgvs.split(';')
                phgvs_list = phgvs.split(';')

                for chgvs in chgvs_list:
                    key = '{gene}_{chgvs}'.format(**locals())
                    driver_info[key] = 1

                for phgvs in phgvs_list:
                    key = '{gene}_{phgvs}'.format(**locals())
                    driver_info[key] = 1
        
        return driver_info

    
    def get_special_driver_info(self):
        '''
        {
            'delins_gene':
                {
                    'gene1': [ex1, ex20],
                    'gene2': ex19,
                }
        }
        '''
        return yaml.load(open(self.driver_yaml))


    def is_driver_according_chgvs(self, **kw):
        '''根据gene和chgvs判断是否为driver基因
        '''
        key = '{gene}_{chgvs}'.format(**kw)
        if kw['driver_info'].get(key):
            return 1
        return 0


    def is_driver_according_phgvs(self, **kw):
        '''根据gene和phgvs判断是否为driver基因
        '''
        key = '{gene}_{phgvs}'.format(**kw)
        if kw['driver_info'].get(key):
            return 1
        return 0


    def is_driver_according_loss(self, **kw):
        '''根据gene和func判断是否为driver基因
        '''
        if not kw['func'].lower() in ('nonsense', 'frameshift'):
            return 0
        else:
            newfunc = 'loss'

        key = '{gene}_{newfunc}'.format(**locals(), **kw)
        if kw['driver_info'].get(key):
            return 1
        return 0

        
    def is_driver_special_gene_delins(self, **kw):
        # gene, phgvs, exon, special_driver_info):

        '''
        1. 部分氨基酸变异类型为delins类型,部分基因状态为driver基因
        传入库文件,不在代码进行判读
        '''
        if not 'delins' in kw['phgvs']:
            return 0

        driver_info = kw['special_driver_info']['delins']
        
        if kw['gene'] in driver_info and kw['exon'] in driver_info[kw['gene']]:
            return 1
        return 0


    def is_driver_special_gene_func(self, **kw):
        #gene, func, exon, special_driver_info):
        '''        
        2. 部分基因->某些func类型的->某些exon， 属于driver基因
        '''
        if not kw['func'] in kw['special_driver_info']['special_func']:
            return 0

        driver_info = kw['special_driver_info'][kw['func']]
        
        if kw['gene'] in driver_info and kw['exon'] in driver_info[kw['gene']]:
            return 1
        return 0


    def is_driver_tp53_special_site(self, **kw):
        '''tp53某个区间内位点属于driver基因
        '''
        if not 'delins' in kw['phgvs'] and not kw['gene'] == 'TP53':
            return 0
        
        if re.search(r'c.(\d+).*', kw['chgvs']):
            chgvs_position = re.search(r'c.(\d+).*', kw['chgvs']).group(1)

            if 102 <= int(chgvs_position) <= 292:
                return 1
        return 0


    def is_driver_gene(self, **kw):
        return any(
            (self.is_driver_according_chgvs(**kw),
            self.is_driver_according_phgvs(**kw),
            self.is_driver_according_loss(**kw),
            self.is_driver_special_gene_delins(**kw),
            self.is_driver_special_gene_func(**kw),
            self.is_driver_tp53_special_site(**kw)
            )
        )


    def start(self):

        driver_info = self.get_driver_info()
        special_driver_info = self.get_special_driver_info()

        with open(self.infile, 'r') as fr, open(self.result, 'w') as fw:
            for line in fr:
                if line.startswith('Scale'): #GG
                    head = line
                    fw.write(line)
                    head_index = utils.get_head_index(head)
                    continue
                linelist = line.strip('\n').split('\t')
                gene = linelist[head_index['hugo_symbol']] 
                chgvs = linelist[head_index['hgvsc']] 
                phgvs = linelist[head_index['hgvsp_short']] 
                func = linelist[head_index['bgi_function']] 
                exon = linelist[head_index['exon']] 
                case_freq = linelist[head_index['case_var_freq']]

                # 利用bgicg结果进行测试
                # gene = linelist[head_index['#gene']] # gene
                # chgvs = linelist[head_index['chgvs']] # cHGVS
                # phgvs = linelist[head_index['phgvs']] # pHGVS
                # func = linelist[head_index['function']] # Function
                # exon = linelist[head_index['exin_id']] # ExIn_ID
                # case_freq = linelist[head_index['case_var_freq']]

                kw = {
                    'gene': gene,
                    'chgvs': chgvs,
                    'phgvs': phgvs,
                    'func': func,
                    'exon': exon,
                    'driver_info': driver_info,
                    'special_driver_info': special_driver_info
                }
                # 进行TMB字段判断
                if float(case_freq) < self.tmb_freq:
                    linelist[head_index['tmb_type']] = 'noTMB'  # 
                elif func == 'span' and re.search(r'-EX1$', exon):
                    linelist[head_index['tmb_type']] = 'noTMB'
                elif self.is_driver_gene(**kw):
                    linelist[head_index['tmb_type']] = 'Driver'
                else:
                    linelist[head_index['tmb_type']] = 'TMB'
 
                fw.write('{}\n'.format('\t'.join(linelist)))




if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='TMB计算代码测试单例')
    parser.add_argument('--infile')
    parser.add_argument('--result')
    parser.add_argument('--driver_list')
    parser.add_argument('--driver_yaml')
    parser.add_argument('--tmb_freq')
    parser.add_argument('--chip_size')

    args = vars(parser.parse_args())

    tmb = TMB(args).start()




