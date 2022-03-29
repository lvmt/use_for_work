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

from simplify import headers
from simplify import  get_vcf_info
from utils import utils
from simplify.transver_function import TransverFunction




class SimpleVep:
    '''
    vep注释程序简化及变异频率添加
    '''
    def parser_vep_simplify(self, parser):
        parser.add_argument(
            '--vep_annotation', help='vep注释原始结果',
        )
        parser.add_argument(
            '--result', help='简化结果文件',
        )
        parser.add_argument(
            '--transcript_data', help='基因-转录本对应文件'
        )
        parser.add_argument(
            '--vep_function_yaml', help='vep function优先级配置文件'
        )
        parser.add_argument(
            '--hg19', help='hg19参考基因组'
        )
        parser.add_argument(
            '--vcf', help='vcf文件,添加变异支持数信息'
        )
        parser.add_argument(
            '--vcftype', help='vcf文件类型'
        )
        parser.set_defaults(func=self.start)
        
    
    def get_tran_relation(self, **args):
        '''
        input:
            self.transcript: gene transcript
        output:
            list: [trans1, trans2]
        '''
        tran_relation = []
        with utils.safe_open(args['transcript_data'], 'r') as fr:
            for line in fr:
                line = line.strip('\n')
                if line.startswith('#'):
                    head_index = utils.get_head_index(line)
                    continue
                linelist = line.strip('\n').split('\t')
                # gene = linelist[head_index['#gene']]
                tran = linelist[head_index['transcript']]
                tran_relation.append('{tran}'.format(**locals()))
        
        return tran_relation


    # def add_freq(self, vcf_info, upload_variation):
    #     '''
    #     如果提供了vcf文件，则给每个变异添加对应的支持数信息
    #     vcf_info: vcf大字典
    #     '''
    #     return vcf_info[upload_variation]


    def start(self, **args):
        '''
        程序运行主函数
        '''
        other_info = {}  # 用于更新其他字段

        if args['vcf']:
            vcf_info = get_vcf_info.HandleVcf(args['vcf'], args['vcftype']).start()
            
        tran_relation = self.get_tran_relation(**args)

        with open(args['vep_annotation'], 'r') as fr, open(args['result'], 'w') as fw:
            fw.write('{}\n'.format('\t'.join(headers.HEAD().update_head().keys())))
            for line in fr:
                if line.startswith('##'):
                    continue
                elif line.startswith('#Uploaded_variation'):
                    head = line
                    head_index = utils.get_head_index(head)
                    continue

                linelist = line.strip('').split('\t')
                gene = linelist[head_index['symbol']] # TERT
                transcript = linelist[head_index['feature']] # NM_198253.3
                
                #提取基因指定的转录本注释信息
                if not '{transcript}'.format(**locals()) in tran_relation:
                    continue

                #获取需要的信息
                row = headers.HEAD()
                ## 可以直接提取的信息
                upload_variation = linelist[head_index['#Uploaded_variation'.lower()]] # chr5_1295229_-/A
                location = linelist[head_index['location']] # chr5:1295187-1295188
                vep_function = linelist[head_index['consequence']] # missense_variant
                strand = linelist[head_index['strand']] # -1
                strand = '+' if strand == '1' else '-'
                protein = linelist[head_index['ensp']] # NP_937983.2
                sift = linelist[head_index['sift']] # tolerated(0.05)
                polyphen = linelist[head_index['polyphen']]
                exon_info = linelist[head_index['exon']] # 2/19 or -
                intro_info = linelist[head_index['intron']]
                chgvs = linelist[head_index['hgvsc']] # NM_198253.3:c.77C>T
                phgvs = linelist[head_index['hgvsp']] # NP_937983.2:p.Thr26Met
                tert = linelist[head_index['tert']] # 只有tert的启动子区域有
                clinvar = linelist[head_index['clinvar_clnsig']]
                rs = linelist[head_index['existing_variation']]
                bl_muttype = linelist[head_index['variant_class']]  
                af = linelist[head_index['af']]
                afr_af = linelist[head_index['afr_af']]
                amr_af = linelist[head_index['amr_af']]
                eas_af = linelist[head_index['eas_af']]
                eur_af = linelist[head_index['eur_af']]
                sas_af = linelist[head_index['sas_af']]
                aa_af = linelist[head_index['aa_af']]
                ea_af = linelist[head_index['eas_af']]
                gnomad_af = linelist[head_index['gnomad_af']]
                gnomad_afr_af = linelist[head_index['gnomad_afr_af']]
                gnomad_amr_af = linelist[head_index['gnomad_amr_af']]
                gnomad_asj_af = linelist[head_index['gnomad_asj_af']]
                gnomad_eas_af = linelist[head_index['gnomad_eas_af']]
                gnomad_fin_af = linelist[head_index['gnomad_fin_af']]
                gnomad_nfe_af = linelist[head_index['gnomad_nfe_af']]
                gnomad_oth_af = linelist[head_index['gnomad_oth_af']]
                gnomad_sas_af = linelist[head_index['gnomad_sas_af']]

                # 需要进行处理获取的信息
                hgvsc = utils.simplify_hgvsc(gene, chgvs, tert)
                hgvsp = utils.simplify_hgvsp(phgvs)  # p.Lys872_Thr874delinsAsnTer
                hgvsp_short = utils.get_oneletter_hgvsp(hgvsp) # p.K872_T874delinsN*
                exon_id = utils.get_exon_id(exon_info, intro_info)
                _chr, start, end = utils.get_chr_start_end_from_location(location)
                ref, alt = utils.get_ref_alt_from_upload_variation(upload_variation)
                muttype = utils.get_muttype(ref, alt)
                genotype = utils.get_genotype(ref, alt, strand)
                flank = utils.get_flank_according_upload_variation(upload_variation, args['hg19'])
                vep_simple_function = TransverFunction(**args).simplify_function(vep_function, tert, gene)
                vep2bgicg_function = TransverFunction(**args).vep2bgi(vep_simple_function, hgvsc, hgvsp, ref, alt, exon_id)

                ## 存在特殊情况,span,跨越整个内含子,但是phgvs还存在注释信息,这种是错误的
                ## 针对这种情况，需要对span类型的phgvs赋空值
                if vep2bgicg_function == 'span' and (not hgvsp == '-'):
                    hgvsp = '-'
                    hgvsp_short = '-'

                # 更新row
                row.gene = gene 
                row.chgvs = hgvsc
                row.phgvs = hgvsp
                row.phgvs_shoft = hgvsp_short
                row.exon_id = exon_id
                # row.tert = tert
                row.vep_function = vep_function
                row.vep_simple_function = vep_simple_function
                row.vep2bgicg_function = vep2bgicg_function
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
                row.af = af
                row.afr_af = afr_af
                row.amr_af = amr_af
                row.eas_af = eas_af
                row.eur_af = eur_af
                row.sas_af = sas_af
                row.aa_af = aa_af
                row.ea_af = ea_af
                row.gnomad_af = gnomad_af
                row.gnomad_afr_af = gnomad_afr_af
                row.gnomad_amr_af = gnomad_amr_af
                row.gnomad_asj_af = gnomad_asj_af
                row.gnomad_eas_af = gnomad_eas_af
                row.gnomad_fin_af = gnomad_fin_af
                row.gnomad_nfe_af = gnomad_nfe_af
                row.gnomad_oth_af = gnomad_oth_af
                row.gnomad_sas_af = gnomad_sas_af

                if args['vcf']:
                    freq_tag = vcf_info[upload_variation]
                    other_info.update(freq_tag)
                info = row.update_head(**other_info)
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






