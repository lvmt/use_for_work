#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   filt_add_cismut.py
@Time    :   2022/01/07 10:28:58
@Author  :   lvmt
@Version :   1.0
'''

# here put the import lib
# 将顺势突变的结果整合至filt中

'''
输入文件：
1. snv/filt_variation
2. 顺势结果文件
3. 顺势raw.vcf
4. T790_T797 check 文件

'''

from os import read
import re
from collections import defaultdict



class GetFinalResult:

    def __init__(self, args):
        self.cismutation = args['cismutation']
        self.cismut_vcf = args['cismut_vcf']
        self.filt = args['filt']
        self.t790 = args['t790']
        self.result = args['result']
        self.header = ''


    def get_head_index(self, head):
        head_info = {}
        if not isinstance(head, list):
            head = head.strip('\n').split('\t')

        for index,item in enumerate(head):
            head_info[item.lower()] = index

        return head_info


    def handle_cismutation(self):
        '''处理顺势突变文件
        info = {'chr_start_ref_alt': linelist}
        '''
        cismut_info = {}
        with open(self.cismutation, 'r') as fr:
            for line in fr:
                if line.startswith('Scale'):
                    head_index = self.get_head_index(line)
                    continue
                linelist = line.strip('\n').split('\t')
                gene = linelist[head_index['hugo_symbol']]
                chgvs = linelist[head_index['hgvsc']]
                _chr = linelist[head_index['chromosome']]
                pos = linelist[head_index['start_position']]
                ref = linelist[head_index['reference_allele']]
                alt = linelist[head_index['allele']]
                if gene == 'TNFRSF14' and chgvs == 'c.-65T>A':
                    continue

                key = '{_chr}-{pos}-{ref}-{alt}'.format(**locals())
                cismut_info[key] = linelist

        return cismut_info


    def handle_cismut_vcf(self):
        '''记录参与顺势合并的原始变异
        '''
        origin_mut_of_merge_cis = []
        with open(self.cismut_vcf, 'r') as fr:
            for line in fr:
                if line.startswith('#'):
                    continue
                linelist = line.strip('\n').split('\t')
                info = linelist[7] # MergeMut=chr7-151935798-G-T:chr7-151935799-A-G;
                if not 'MergeMut' in info:
                    continue
                muta, mutb = re.search(r'MergeMut=([a-zA-Z0-9\-]+):([a-zA-Z0-9\-]+);', info).groups()
                origin_mut_of_merge_cis.append(muta)
                origin_mut_of_merge_cis.append(mutb)
        return origin_mut_of_merge_cis


    def handle_filt(self, origin_mut_of_merge_cis):
        filt_info = {}
        with open(self.filt, 'r') as fr:
            for line in fr:
                if line.startswith('Scale'):
                    head_index = self.get_head_index(line)
                    self.header = line.strip('\n').split('\t')
                    continue
                linelist = line.strip('\n').split('\t')
                _chr = linelist[head_index['chromosome']]
                pos = linelist[head_index['start_position']]
                ref = linelist[head_index['reference_allele']]
                alt = linelist[head_index['allele']]
                key = '{_chr}-{pos}-{ref}-{alt}'.format(**locals())

                if key in origin_mut_of_merge_cis:  # 丢弃参与顺势合并的原始变异位点
                    continue

                filt_info[key] = linelist

        return filt_info


    def handle_t790(self, cismut_info, filt_info):
        t790_info = defaultdict(str)
        t790_totalread_info = defaultdict(int)

        gene_chgvs_list = []
        for d in (cismut_info, filt_info):
            for key in d:
                linelist = d[key]
                gene = linelist[1]
                chgvs = linelist[2]
                gene_chgvs_list.append('{gene}:{chgvs}'.format(**locals()))

        with open(self.t790, 'r') as fr:
            for line in fr:
                if 'Undetect' in line:
                    continue
                tlinelist = line.strip('\n').split('\t')
                muta = tlinelist[0] # EGFR:c.2369C>T
                mutb = tlinelist[1] # EGFR:c.2389T>A
                ty = tlinelist[2]   # Trans{c.2369C>T}{c.2389T>A}
                support_read = tlinelist[3] # 7|1
                if not (muta in gene_chgvs_list and mutb in gene_chgvs_list):
                    # 只有muta 和 mutb同时存在于结果文件中,才进行后续分析
                    continue

                if ty.startswith('Cis'):
                    t790_info[muta] += '{ty}:{support_read};'.format(**locals())
                    t790_info[mutb] += '{ty}:{support_read};'.format(**locals())
                    t790_totalread_info[muta] += int('{support_read}'.format(**locals()))
                    t790_totalread_info[mutb] += int('{support_read}'.format(**locals()))

                elif ty.startswith('Trans'):
                    read_list = list(map(int, support_read.split('|')))
                    if read_list[0] > 5:
                        t790_info[muta] += '{ty}:{read_list[0]};'.format(**locals())
                    if read_list[1] > 5:
                        t790_info[mutb] += '{ty}:{read_list[1]};'.format(**locals())
                    t790_totalread_info[muta] += read_list[0]
                    t790_totalread_info[mutb] += read_list[1]

        return t790_info, t790_totalread_info


    def start(self):
        cismut_info = self.handle_cismutation()
        origin_mut_of_merge_cis = self.handle_cismut_vcf()
        filt_info = self.handle_filt(origin_mut_of_merge_cis)
        t790_info, t790_totalread_info = self.handle_t790(cismut_info, filt_info)
        head_index = self.get_head_index(self.header)

        with open(self.result, 'w') as fw:
            fw.write('{}\n'.format('\t'.join(self.header)))
            for info in (cismut_info, filt_info):
                for key in info:
                    linelist = info[key]
                    gene = linelist[1]
                    chgvs = linelist[2]
                    var_freq = linelist[head_index['case_var_freq']]
                    t790_key = '{gene}:{chgvs}'.format(**locals())
                    if t790_key in t790_info:  #  'ALK:c.3604G>A' => 'Trans{c.3604G>A}{c.3617C>T}:6;Cis{c.3586C>A}{c.3604G>A}:4;'
                        tag_info = ''
                        value_list = t790_info[t790_key].strip(';').split(';')
                        for item in value_list:
                            tag, read = item.split(':')
                            freq = float(var_freq) * float(read) / float(t790_totalread_info[t790_key])
                            tag_info += '{tag}:{freq};'.format(**locals())

                        linelist[head_index['bgi_newfunction']] = tag_info

                    fw.write('{}\n'.format('\t'.join(linelist)))



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='add cismutation to filt')
    parser.add_argument('--filt', help='filt文件,1.snv/filt_variation')
    parser.add_argument('--cismutation', help='顺势变异的结果文件')
    parser.add_argument('--cismut_vcf', help='顺势结果的原始vcf文件')
    parser.add_argument('--t790', help='T790_check文件')
    parser.add_argument('--result', help='最终的结果文件')

    args = vars(parser.parse_args())

    GetFinalResult(args).start()








