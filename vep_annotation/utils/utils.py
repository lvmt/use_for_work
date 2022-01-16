#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   utils.py
@Time    :   2021/11/12 10:41:23
@Author  :   lvmt 
@Version :   1.0
'''



# here put the import lib
# 数据处理集合

import re 
import subprocess


acid_info = {
    'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D',
    'Cys': 'C', 'Gln': 'Q', 'Glu': 'E', 'Gly': 'G',
    'His': 'H', 'Ile': 'I', 'Leu': 'L', 'Lys': 'K',
    'Met': 'M', 'Phe': 'F', 'Pro': 'P', 'Ser': 'S',
    'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V',
    'Ter': '*'
}


#######################################
############ 文件处理函数 ##############
def safe_open(filename, mode='r'):
    if filename.endswith('.gz'):
        import gzip
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)


def get_head_index(headlist):
    head_index = {}
    if not isinstance(headlist, list):
        headlist = headlist.strip('\n').split('\t')

    for index,item in enumerate(headlist):
        head_index[item.lower()] = index

    return head_index


def reverse_complement(base):
    base = base[::-1].upper()
    return base.replace('A', 't').\
                replace('C', 'g').\
                replace('T', 'a').\
                replace('G', 'c').upper()



########################################
##########字段处理函数###################
def get_base(chr, start, end, hg19):
    cmd = 'samtools faidx {hg19} {chr}:{start}-{end}'.format(**locals())
    result = subprocess.getoutput(cmd).split('\n')[-1]
    return result


def get_flank_according_upload_variation(upload_variation, hg19):
    '''获取变异的侧翼序列
    pos：为变异发生的真实位置信息1base
    SNV: POS-2,POS-1 | POS+1,POS+2
    INS: POS-2,POS-1 | POS, POS+1
    DEL: POS-2,POS-1 | (POS + len -1 ) + 1, (POS + len -1 ) + 2
    INDEL: POS-2, POS-1 | (POS + len -1 ) + 1, (POS + len -1 ) + 2
    '''

    itemlist = re.split('_|/', upload_variation)
    _chr = itemlist[0]
    pos = int(itemlist[1])
    ref = itemlist[2]
    alt = itemlist[3]

    muttype = get_muttype(ref, alt)

    lstart = pos - 2
    lend = pos - 1

    if muttype == 'snv':
        rstart = pos + 1
        rend = pos + 2
    elif muttype == 'ins':
        rstart = pos
        rend = pos + 1
    else:
        rstart = (pos + len(ref) - 1) + 1
        rend = (pos + len(ref) - 1) + 2

    lbase = get_base(_chr, lstart, lend, hg19)
    rbase = get_base(_chr, rstart, rend, hg19)
    return '.'.join((lbase, rbase))


def get_muttype(ref, alt):
    if ref == '-':
        return 'ins'
    elif alt == '-':
        return 'del'
    elif len(ref) + len(alt) == 2:
        return 'snv'
    else:
        return 'delins'


def get_muttype2(ref, alt):
    '''根据vcf文件的ref和alt判断变异类型
    '''
    if ref[0] != alt[0] and len(ref) + len(alt) > 2:
        return 'delins'
    elif len(ref) == 1 and len(alt) > 1:
        return 'ins'
    elif len(alt) == 1 and len(ref) > 1:
        return 'del'
    else:
        return 'snv'


def modify_pos_ref_alt(pos, ref, alt):
    '''根据indel vcf的pos,ref,alt，生成新的pos,ref,alt
    更新后的pos,ref,alt和vep的Uploaded_variation中信息一致
    del: pos + 1
    ins: pos + 1
    delins: pos 
    '''
    muttype = get_muttype2(ref, alt)
    if muttype == 'del':
        pos = int(pos) + 1
        ref = ref[1:]
        alt = '-'
    elif muttype == 'ins':
        pos = int(pos) + 1
        ref = '-'
        alt = alt[1:]

    return pos, ref, alt       


def get_chr_start_end_from_location(location):

    '''根据vep的location字段获取变异的chr,start,end

    snv: chr5_1295187_G/A	chr5:1295187
    ins: chr5_1295188_-/GGG	chr5:1295187-1295188
    del: chr5_1295233_CTGG/-	   chr5:1295233-1295236
    delins: chr5_1295232_G/TT	   chr5:1295232
    delins: chr5_1295232_GC/TT   chr5:1295232-1295233
    '''
    linelist = re.split(':|-', location)
    _chr = linelist[0]
    if len(linelist) == 2:
        start = linelist[1]
        end = linelist[1]
    else:
        start = linelist[1]
        end = linelist[2]
    
    return _chr, start, end


def simplify_hgvsc(gene, hgvsc, tert):
    '''处理hgvsc,获取解读需要格式
    争对tert启动子区域需要特殊处理
    input: NM_198253.3:c.77C>T
    output: c.77C>T
    tert: TERT字段
    '''
    if gene == 'TERT' and tert.startswith('NM'):
        return tert.split(':')[-1]
    return hgvsc.split(':')[-1]


def simplify_hgvsp(hgvsp):
    '''处理hgvsp
    input: NP_937983.2:p.Thr26del
    output：p.Thr26del
    '''
    return hgvsp.split(':')[-1]


def get_oneletter_hgvsp(hgvsp):
    '''获取单字母的hgvsp结果
    '''
    for acid in acid_info:
        if acid in hgvsp:
            hgvsp = hgvsp.replace(acid, acid_info[acid])
    return hgvsp


def get_exon_id(exon_info, intro_info):
    '''格式处理
    根据exon和intro列信息进行判断
    input: 2/16
    特殊：对于span位点，会同时存在exon和intro：2021-1201
    特殊：存在大片段跨区变异：2021-1213
    15/18 14-15/17
    '''
    exon_info = exon_info.split('/')[0]
    intro_info = intro_info.split('/')[0]
    
    if exon_info == intro_info == '-':
        return '-'
    elif not exon_info == '-' and not intro_info == '-':
        return 'IVS{intro_info}-EX{exon_info}'.format(**locals())
    elif not exon_info == '-':
        return 'EX{exon_info}'.format(**locals())
    elif not intro_info == '-':
        return 'IVS{intro_info}'.format(**locals())


def get_ref_alt_from_upload_variation(upload_variation):
    '''利用upload_variation字段,获取变异的ref和alt
    '''
    itemlist = upload_variation.split('_')[-1].split('/')
    ref = itemlist[0]
    alt = itemlist[1]
    return ref, alt


def get_genotype(ref, alt, strand):
    '''
    如果是负链基因，需要做对应的反向互补操作
    snv: A/T
    ins: ./CCGG
    del: GAA/.
    delins: AG/TA
    2021-11-20好像不需要反向互补
    '''
    # if int(strand) == -1:
    #     ref = reverse_complement(ref)
    #     alt = reverse_complement(alt)
    return '{ref}/{alt}'.format(**locals())


def get_gtype(GT):
    '''根据vcf的GT字段获取，分为杂合和纯合
    0/1: 杂合
    1/1：纯合
    '''
    if GT == '0/1':
        return 'het'
    elif GT == '1/1':
        return 'hom'
    else:
        # 可能存在其他类型，开发此函数时，暂未发现
        return 'check!!-{GT}'.format(**locals())


def get_pos_ref_alt_from_vcf(pos, ref, alt):
    '''
    pos：对应maf中的start
    根据vcf格式文件，更新坐标体系及ref/alt
    snv: 不变
    indel: 123 CT C -> 124 T -
    ins: 123 C CTT -> 123 - TT
    delins: 123 C TT -> 123 C TT
    delins: 123 CA TT -> 123 CA TT
    '''
    muttype = get_muttype2(ref, alt)
    if muttype in ('snv', 'delins'):
        return pos, ref, alt
    elif muttype == 'del':
        return int(pos) + 1, ref[1:], '-'
    elif muttype == 'ins':
        return pos, '-', alt[1:]

