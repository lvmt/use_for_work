#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   combine_hot_filt.py
@Time    :   2021/12/06 15:49:48
@Author  :   lvmt 
@Version :   1.0
'''

# here put the import lib
'''
数据合并
目的：将hot中的结果，添加到filt中
hot位点全都是snv
处理：
    1. hot中需要剔除filt中存在的结果
    2. hot中需要剔除遗传性位点中存在的结果
'''


import pandas as pd
from pandas.core.reshape.merge import merge 



class CombineFiltHot:

    def __init__(self, args):
        self.filt = args['filt']
        self.hot = args['hot']
        self.here = args['here_tumor']
        self.result = args['result']


    def handle_filt(self):
        filt_df = pd.read_csv(self.filt, sep='\t', dtype=str)
        filt_df['cpra'] = filt_df['Chromosome'] + '-' + filt_df['End_Position'] + '-' + filt_df['Reference_Allele'] + '-' + filt_df['Allele']
        return filt_df


    def handle_hot(self):
        '''注意坐标体系
        chr1	115252202	115252203	G	T
        '''
        hot_df = pd.read_csv(self.hot, sep='\t', dtype=str)
        hot_df.rename(
            columns = {
                '#Case_score': 'Scale',
                'Ctrl_score': 'Scale_ctrl',
                'Case_ratio': 'Case_var_freq',
                'Chr': 'Chromosome',
                'Start': 'Start_Position',
                'Stop': 'End_Position',
                'Ref': 'Reference_Allele',
                'Call': 'Allele',
                'Gene Symbol': 'Hugo_Symbol',
                'ExIn_ID': 'EXON',
                'Function': 'BGI_Function',
                'Transcript': 'Transcript_ID',
                'Protein': 'ENSP',
                'cHGVS': 'HGVSc',
                'pHGVS': 'HGVSp_Short',
                'Uniport_Position(s)': 'BGI_Uniport_Position(s)',
                'Uniport_feature_key': 'BGI_Uniport_feature_key',
                'Uniport_description': 'BGI_Uniport_description',
                'cosmic': 'Cosmic',
               }, inplace=True
        )
        
        # 修改坐标体系
        hot_df['Start_Position'] = hot_df['End_Position']

        hot_df['cpra'] = hot_df['Chromosome'] + '-' + hot_df['End_Position'] + '-' + hot_df['Reference_Allele'] + '-' + hot_df['Allele']

        hot_df['Case_pos_dep'] = hot_df['Case_reads'].map(lambda x: x.split('/')[1]) # 9/2799
        hot_df['Case_var_readsNum'] = hot_df['Case_reads'].map(lambda x: x.split('/')[0])
        hot_df['Case_ref_readsNum'] = [str(int(x) - int(y)) if x != '*' else '*' for x, y in zip(hot_df['Case_pos_dep'], hot_df['Case_var_readsNum'])]

        hot_df['Ctrl_pos_dep'] = hot_df['Ctrl_reads'].map(lambda x: x.split('/')[1] if '/' in x else x)
        hot_df['Ctrl_var_readsNum'] = hot_df['Ctrl_reads'].map(lambda x: x.split('/')[0] if '/' in x else x)
        hot_df['Ctrl_ref_readsNum'] = [str(int(x) - int(y)) if x != '*' else '*' for x, y in zip(hot_df['Ctrl_pos_dep'], hot_df['Ctrl_var_readsNum'])]

        return hot_df


    def handle_here(self):
        '''注意坐标体系
        snv chr17   29553484       29553484       G       T | chr17-29553484-C-T (cpra)
        '''
        here_df = pd.read_csv(self.here, sep='\t', dtype=str)
        return here_df


    def get_target_hot(self, hot_df, filt_df, here_df):
        
        hot_df = hot_df[~hot_df['cpra'].isin(filt_df['cpra'].tolist())]  # 剔除filt中和hot相同的cpra位点
        hot_df = hot_df[~hot_df['cpra'].isin(here_df['cpra'].tolist())]  # 剔除here中和hot相同的cpra位点

        hot_df['HGVSp_Short'] = hot_df['HGVSp_Short'].map(lambda x: x.split('|')[0])
        # hot_df['score'] = hot_df['Scale'].map(lambda x: x.split('/')[0])

        # 肺癌验证位点
        lung_hot_judge = (hot_df['Case_var_freq'].astype(float) < 0.6) & (hot_df['Hugo_Symbol'] == 'EGFR') & (hot_df['HGVSp_Short'].str.contains(r'L858R|T790M'))
        hot_df['Scale'][lung_hot_judge] = hot_df['Scale'][lung_hot_judge].map(lambda x: '肺癌验证.' + x) ##只包含snv，是否需要19del
        hot_df['Scale'] = 'Hot.' + hot_df['Scale']

        # 增加字段
        hot_df['Codons'] = hot_df['Reference_Allele'] + '/' + hot_df['Allele']
        hot_df['Variant_Type'] = 'SNV'

        return hot_df


    def start(self):
        hot_df = self.handle_hot()
        filt_df = self.handle_filt()
        here_df = self.handle_here()

        target_hot = self.get_target_hot(hot_df, filt_df, here_df)
        del filt_df['cpra']

        merge_df = pd.merge(filt_df, target_hot, how='outer')[filt_df.columns].fillna('*')
        merge_df.to_csv(self.result, index=None, sep='\t')

        # print(merge_df)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='将hot变异添加至filt中')
    parser.add_argument('--filt', help='filt文件')
    parser.add_argument('--hot', help='hot文件')
    parser.add_argument('--here_tumor', help='遗传性位点结果')
    parser.add_argument('--result', help='最终的结果文件')

    args = vars(parser.parse_args())

    CombineFiltHot(args).start()


