#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   transver_function.py
@Time    :   2021/11/15 14:45:09
@Author  :   lvmt 
@Version :   1.0
'''

# here put the import lib
'''将vep原始注释功能转换为对应的BGICG功能
1. 先根据规则注释功能进行转换
2. 在根据功能yaml文件,进行替换
'''


'''
注释结果中2类特殊function需要优先处理
splice_region_variant: 遇到这个func,直接从func_list中del
protein_altering_variant：第一优先级,需要根据氨基酸的phgvs进行进一步判断
'''


from ast import arg
from distutils.spawn import spawn
from sqlalchemy import desc
import yaml
import re
# from simplify import simplify_vep_annotation
# import utils  




class TransverFunction:

    def __init__(self, **args):
        self.vep_function_yaml = args['vep_function_yaml']


    def simplify_function(self, vep_function, tert, gene):
        '''
        优先处理tert 启动子区域问题
        vep的注释结果可能存在多个,需要简化成为一个
        '''
        #因为splice_region_variant是个陪衬,先清除掉
        vep_function = vep_function.replace('splice_region_variant', '').strip(',')
        func_list = list(filter(str, vep_function.split(','))) # 去除多余的分隔符 &

        if len(func_list) == 1 and gene == 'TERT' and tert.startswith('NM'):
            return 'promter'
        elif len(func_list) == 1:
            return vep_function
        elif len(func_list) > 1:  # 及时func_list=3, 也不一定是span
            return self.get_the_better_func(func_list)

    
    def get_the_better_func(self, func_list):
        '''
        1.有时候2个注释也会是span，因为跨过了不同的编码区域
        2.根据function的优先级,返回优先级较高的func信息
        '''
        region_set = set()  # 存储func的exon or intro信息
        func_region_info = yaml.load(open(self.vep_function_yaml))['FuncRegion']
        func_priority = yaml.load(open(self.vep_function_yaml))['Priority']

        for func in func_list:
            region_set.add(func_region_info[func])

        if len(region_set) > 1:
            return 'span'
        else:
            # 可能存在3个或者多个相同的region
            # frameshift/start-lost/start-retained
            # 选取其中优先级数值最小的func作为最优func
            better_func = func_list[0]
            for func in func_list:
                if func_priority[func] < func_priority[better_func]:
                    better_func = func
            
            return better_func


    def vep2bgi(self, simplify_vep_func, hgvsc, hgvsp, ref, alt):
        '''将vep的func转化为bgi的func
        优先处理protein_altering_variant
        优先处理coding_sequence_variant
        对于span需要进行再次矫正
        对于splice,由于存在特殊情况,因此也需要进行再次矫正
        '''
        vep2bgi_info = yaml.load(open(self.vep_function_yaml))['vep2BGI']

        if simplify_vep_func == 'protein_altering_variant':
            return self.handle_protein_altering_variant_from_hgvsp(hgvsp, ref, alt)
        elif simplify_vep_func == 'coding_sequence_variant':
            # return 'coding_sequence_variant'
            return self.handle_coding_sequence_variant_from_chgvs_phgvs(hgvsc, hgvsp, ref, alt)
        elif simplify_vep_func == 'span':
            return self.correct_span(hgvsp, ref, alt)
        elif simplify_vep_func.startswith('splice'):
            return self.correct_splice(hgvsc)
        else:
            return vep2bgi_info[simplify_vep_func]


    def handle_protein_altering_variant_from_hgvsp(self, hgvsp, ref, alt):
        '''根据phgvs判断bgi的function结果
        现根据phgvs判断:stop_gain,frameshift
        在根据ref，alt判断: 
        cds-ins, cds-del,cds-indel  
        '''
        if 'fs' in hgvsp:
            return 'frameshift'
        elif hgvsp.endswith('Ter'):
            # stop_gained
            return 'nonsense'
        elif 'delins' in hgvsp:  #还需进一步判断,解读逻辑白帅帅
            if len(ref) > len(alt):
                return 'cds-del'
            elif len(ref) < len(alt):
                return 'cds-ins'
            else:
                return 'missense'  # 氨基酸已经发生变化，20220217
        elif 'del' in hgvsp:
            return 'cds-del'
        elif 'ins' in hgvsp:
            return 'cds-ins'
        elif '=' in hgvsp:
            return 'coding-synon'
        elif re.search(r'[A-Za-z]+\d*[A-Za-z]+', hgvsp):
            return 'missense'
            # chr8_145738769_G/A
            # splice_acceptor_variant,missense_variant
        else:
            return 'null'


    def handle_coding_sequence_variant_from_chgvs_phgvs(self, hgvsc, hgvsp, ref, alt):
        '''
        coding_sequence_variant
        决策就是根据chgvs和phgvs进行判断
        若存在phgvs,则利用handle_proten***这个函数进行判断
        否则按照chgvs进行判断
        c.3727_3727+1delinsAC 
        好像还有更复杂的格式问题，需要解决：20211117
        这个也有部分问题待解决
        ''' 
        if not hgvsp == '-':
            return self.handle_protein_altering_variant_from_hgvsp(hgvsp, ref, alt)
        elif not hgvsc == '-':
            # 因为注释为coding_sequence,
            # 基于chgvs进行判断，可能存在的情况: span, 
            # 内含子区域或者剪切区
            # 更具chgvs进行判断：优先根据_进行判断,是否存在多碱基的del等
            region_tag = set()
            hgvsc = hgvsc.split('c.')[-1].split('_')
            # 只有1个
            if len(hgvsc) == 1:
                # if re.search(r'(\d+)[-][12][^0-9]', hgvsc[0]):
                if re.search(r'(\d+)[-][12]($|[a-zA-Z])', hgvsc[0]):
                    return 'splice-3'
                elif re.search(r'(\d+)[+][12]($|[a-zA-Z])', hgvsc[0]):
                    return 'splice-5'
                elif re.search(r'(\d+)[+-][3-9]*', hgvsc[0]):
                    return 'intro'
                else:
                    return 'null'
            else:
                # 特例：NM_001081640.2:c.3727+1_3728-1insT
                for item in hgvsc:
                    if '-' in item or '+' in item:
                        region_tag.add('intro')
                    else: 
                        region_tag.add('exon')

                if len(region_tag) > 1:
                    return 'span'
                else:
                    # 多个区间信息,没有氨基酸注释,不会是fs,missense
                    if re.findall('(\d+)[+-][12]($|[a-zA-Z])', hgvsc[0]):
                        return 'splice-5'
                    else:
                        return 'intro'
            

            
    def correct_span(self, hgvsp, ref, alt):
        '''
        vep中出现一类特殊案例
        splice_acceptor_variant&frameshift_variant&intron_variant
        基于上述的逻辑判断，最终的功能会注释为span，但是在回查注释结果时，phgvs是由结果的
        因此对于注释为span的类型，我们需要进行再次矫正
        '''
        if hgvsp == '-':
            return 'span'
        else:
            return self.handle_protein_altering_variant_from_hgvsp(hgvsp, ref, alt)


    def correct_splice(self, hgvsc):
        # 大数据量测试过程中发现, 存在一下情况
        # splice_donor_variant,intron_variant -> NM_006437.3:c.3285+3_3285+5del
        # NM_014727.1:c.3059-11del（不可以匹配成-1）！！！
        if re.findall(r'[+][1-2][^0-9]', hgvsc):
            return 'splice-5'
        elif re.findall(r'[-][1-2][^0-9]', hgvsc):
            return 'splice-3'
        else:
            return 'intron'



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='test')
    parser.add_argument('--infile')
    parser.add_argument('--result')
    parser.add_argument('--vep_function_yaml', help='func 配置文件')
    args = vars(parser.parse_args())

    infile = args['infile']
    result = args['result']

    with open(infile, 'r') as fr, open(result, 'w') as fw:
        for line in fr:
            linelist = line.strip('\n').split('\t')
            if line.startswith('#'):
                head_index = utils.get_head_index(linelist)
                fw.write('{}\tvep_simple\tbgi_func\n'.format('\t'.join(linelist)))
                continue

            nm = linelist[head_index['feature']]
            if not nm.startswith('NM'):
                continue

            gene = linelist[head_index['symbol']]
            upload_variation = linelist[head_index['#Uploaded_variation'.lower()]]
            vep_function = linelist[head_index['consequence']]
            chgvs = linelist[head_index['hgvsc']]
            phgvs = linelist[head_index['hgvsp']]
            tert = linelist[head_index['tert']]
            

            hgvsc = utils.simplify_hgvsc(gene, chgvs, tert)
            hgvsp = utils.simplify_hgvsp(phgvs) 
            ref, alt = utils.get_ref_alt_from_upload_variation(upload_variation)
            vep_simple_function = TransverFunction(**args).simplify_function(vep_function, tert, gene)
            vep2bgicg_function = TransverFunction(**args).vep2bgi(vep_simple_function, hgvsc, hgvsp, ref, alt)

            fw.write('{}\t{}\t{}\n'.format('\t'.join(linelist), vep_simple_function, vep2bgicg_function))




            

