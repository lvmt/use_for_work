#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   get_vep_annotation.py
@Time    :   2021/11/11 14:27:28
@Author  :   lvmt 
@Version :   1.0
'''

# here put the import lib
'''
从vep注释结果中，提取必须的列
'''



import re
import subprocess

from utils import utils
import headers
import get_vcf_info




class SimpleVep:
    '''
    vep注释程序简化及变异频率添加
    '''
    def __init__(self, args):
        self.vep = args['vep_annotation']
        self.transript_database = args['transript_database'] 
        self.vcf = args['vcf']
        self.vcftype = args['vcftype']
        self.hg19 = args['hg19']

    
    def get_tran_relation(self):
        '''
        input:
            self.transcript: gene transcript
        output:
            list: [gene=trans, gene=trans]
        '''
        tran_relation = []
        with utils.safe_open(self.transript_database, 'r') as fr:
            for line in fr:
                if line.startswith('#'):
                    head_index = utils.get_head_index(line)
                    continue
                linelist = line.strip('').split('\t')
                gene = linelist[head_index['#gene']]
                tran = linelist[head_index['transcript']]
                tran_relation.append('{gene}={tran}'.format(**locals()))
        
        return tran_relation


    # def add_freq(self, vcf_info, upload_variation):
    #     '''
    #     如果提供了vcf文件，则给每个变异添加对应的支持数信息
    #     vcf_info: vcf大字典
    #     '''
    #     return vcf_info[upload_variation]


    def start(self):
        '''
        程序运行主函数
        '''
        if self.vcf:
            vcf_info = get_vcf_info.HandleVcf(self.vcf, self.vcftype).start()

        with open(self.vep, 'r') as fr, open('test', 'w') as fw:
            for line in fr:
                if line.startswith('##'):
                    continue
                elif line.startswith('#Uploaded_variation'):
                    head = line
                    head_index = utils.get_head_index(head)
                    continue

                linelist = line.strip('').split('\t')
                #获取需要的信息
                row = headers.HEAD()
                ## 可以直接提取的信息
                upload_variation = linelist[head_index['#Uploaded_variation'.lower()]] # chr5_1295229_-/A
                location = linelist[head_index['location']] # chr5:1295187-1295188
                transcript = linelist[head_index['feature']] # NM_198253.3
                function = linelist[head_index['consequence']] # missense_variant
                strand = linelist[head_index['strand']] # -1
                gene = linelist[head_index['symbol']] # TERT
                protein = linelist[head_index['ensp']] # NP_937983.2
                sift = linelist[head_index['sift']] # tolerated(0.05)
                polyphen = linelist[head_index['polyphen']]
                exon_id = linelist[head_index['exon']] # 2/19 or -
                chgvs = linelist[head_index['hgvsc']] # NM_198253.3:c.77C>T
                phgvs = linelist[head_index['hgvsp']] # NP_937983.2:p.Thr26Met
                tert = linelist[head_index['tert']] # 只有tert的启动子区域有
                clinvar = linelist[head_index['clinvar_clnsig']]
                rs = linelist[head_index['existing_variation']]

                # 需要进行处理获取的信息
                hgvsc = utils.simplify_hgvsc(chgvs)
                hgvsp2 = utils.simplify_hgvsp(phgvs)
                hgvsp = utils.get_oneletter_hgvsp(hgvsp2)
                exon_id = utils.get_exon_id(exon_id)
                _chr, start, end = utils.get_chr_start_end_from_location(location)
                ref, alt = utils.get_ref_alt_from_upload_variation(upload_variation)
                muttype = utils.get_muttype(ref, alt)
                genotype = utils.get_genotype(ref, alt, strand)
                flank = utils.get_flank_according_upload_variation(upload_variation, self.hg19)
                bl_muttype = utils.get_bl_muttype()

                # 更新row
                row.gene = gene 
                row.chgvs = hgvsc
                row.phgvs = hgvsp
                row.phgvs2 = hgvsp2
                row.exon_id = exon_id
                row.vep_function = function
                row.sift = sift
                row.polyphen2 = polyphen
                row.chr = _chr
                row.start = start
                row.end = end 
                row.ref = ref
                row.alt = alt
                row.muttype = muttype
                row.genotype = genotype
                row.transcript = transcript
                row.protein = protein
                row.strand = strand
                row.flank = flank
                row.rs = rs 
                row.bl_muttype = bl_muttype
                row.clinvar = clinvar

                if self.vcf:
                    freq_tag = vcf_info[upload_variation]
                    # print(freq_tag)
                info = row.update_head(**freq_tag)
                fw.write('\t'.join(info.keys()) + '\n')
                fw.write('\t'.join(map(str, info.values())) + '\n')
                # print(info)



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='VEP注释结果简化')
    parser.add_argument('--vep_annotation', help='vep注释原始文件')
    parser.add_argument('--result', help='结果文件')
    parser.add_argument('--transript_database', help='基因对应的转录本文件', required=True)
    parser.add_argument('--hg19', help='计算侧翼序列时需要')
    parser.add_argument('--vcf', help='该文件存在,可以对变异增加变异支持数')
    parser.add_argument('--vcftype', help='vcf文件类型', choices=['snv', 'indel'])

    args = vars(parser.parse_args())

    simple_vep = SimpleVep(args).start()






