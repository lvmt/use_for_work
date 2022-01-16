#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   headers.py
@Time    :   2021/11/11 18:02:30
@Author  :   lvmt 
@Version :   1.0
'''

# here put the import lib
# 存放bgicg默认需要的字段信息



from collections import defaultdict



class HEAD:
    def __init__(self):
        ############################## 
        ######## 基本信息 ############
        self.gene = '*'
        self.chgvs = '*'
        self.phgvs = '*'
        self.phgvs2 = '*'
        self.exon_id = '*'
        self.vep_function = '*'
        self.vep_simple_function = '*'
        self.vep2bgicg_function = '*'
        self.newfunction = '*'
        self.sift = '*'
        self.polyphen2 = '*'
        self.chr = '*'
        self.start = '*'
        self.end = '*'
        self.ref = '*'
        self.alt = '*'
        self.muttype = '*'
        self.genotype = '*'
        self.transcript = '*'
        self.protein = '*'
        self.strand = '*'
        self.flank = '*'
        self.rs = '*'
        self.bl_muttype = '*'
        self.clinvar = '*'
        # 变异频率
        self.af = '*'
        self.afr_af = '*'
        self.amr_af = '*'
        self.eas_af = '*'
        self.eur_af = '*'
        self.sas_af = '*'
        self.aa_af = '*'
        self.ea_af = '*'
        self.gnomad_af = '*'
        self.gnomad_afr_af = '*'
        self.gnomad_amr_af = '*'
        self.gnomad_asj_af = '*'
        self.gnomad_eas_af = '*'
        self.gnomad_fin_af = '*'
        self.gnomad_nfe_af = '*'
        self.gnomad_oth_af = '*'
        self.gnomad_sas_af = '*'
        #######################################
        ############# vcf信息##################
        self.case_depth = '*'
        self.control_depth = '*'
        self.case_var_freq = '*'
        self.case_ref_read = '*'
        self.case_alt_read = '*'
        self.case_alt_positive_read = '*'
        self.case_alt_negative_read = '*'
        self.control_ref_read = '*'
        self.control_alt_read = '*'
        self.control_alt_positive_read = '*'
        self.control_alt_negative_read = '*'
        self.Gtype = '*'
        #####################################
        # uniport 库文件
        self.uniport_position = '*'
        self.uniport_featurn = '*'
        self.uniport_description = '*'
        ######################################
        # cosmic 库文件
        self.cosmic = '*'
        self.cosmic2 = '*'
        ######################################
        # maploc 库文件 变异区带信息
        self.maploc = '*'
        ######################################
        # 变异cds区域 库文件
        self.funcregion = '*'
        ######################################
        # 是否是靶点基因
        self.target_gene = '*'
        ######################################
        # 打分结果文件
        self.case_score = '*'
        self.control_score = '*'
        ######################################
        # TMB 结果
        self.tmb = '*'
        ######################################
        ###########逻辑判断####################
        self.end_exon_check = '*'
        self.gene_extend_check = '*'
        self.splice_affect_check = '*'
        ######################################
        ##############单组织##################
        self.filter_001 = '*'
        self.filter_005 = '*'
        self.interpret = '*'
        # 单组织 天华
        self.in_tianhua = '*'
        self.tianhua_time = '*'
        self.same_position_variant = '*'
        self.tequila_literature = '*'
        self.tequila_class = '*'
        self.tequila_acmg = '*'
        self.primer_design = '*'
     

    def update_head(self, **tags):
        #各个字段默认为*，如果有修改,进行更新
        #tags:其他添加字段
        info = defaultdict()
        info = {
            'Scale': self.case_score,
            'Gene': self.gene,
            'cHGVS': self.chgvs,
            'pHGVS': self.phgvs,
            'pHGVS2': self.phgvs2,
            'EXIn_ID': self.exon_id,
            'Case_var_freq': self.case_var_freq,
            'VEP_Function': self.vep_function,
            'VEP_simple_Function': self.vep_simple_function,
            'vep2bgicg_function': self.vep2bgicg_function,
            'NewFunction': self.newfunction,
            'EndExonCheck': self.end_exon_check,
            'GeneExtendCheck': self.gene_extend_check,
            'Splice_Affect_Check': self.splice_affect_check,
            'Uniport_Position(s)': self.uniport_position,
            'Uniport_feature_key': self.uniport_featurn,
            'Uniport_description': self.uniport_description,
            'cosmic': self.cosmic,
            'Cosmic2': self.cosmic2,
            'SIFT': self.sift,
            'PolyPhen2': self.polyphen2,
            'Case_ref_readsNum': self.case_ref_read,
            'Case_var_readsNum': self.case_alt_read,
            'Case_var_Positive_readsNum': self.case_alt_positive_read,
            'Case_var_Negative_readsNum': self.case_alt_negative_read,
            'Ctrl_ref_readsNum': self.control_ref_read,
            'Ctrl_var_readsNum': self.control_alt_read,
            'Ctrl_var_Positive_readsNum': self.control_alt_positive_read,
            'Ctrl_var_Negative_readsNum': self.control_alt_negative_read,
            'Chr': self.chr,
            'Start': self.start,
            'End': self.end,
            'Ref': self.ref,
            'Alt': self.alt,
            'MapLoc': self.maploc,
            'MutType': self.muttype,
            'Gtype': self.Gtype,
            'Genotype': self.genotype,
            'Transcript': self.transcript,
            'Protein': self.protein,
            'Strand': self.strand,
            'FuncRegion': self.funcregion,
            'Flank': self.flank,
            'Case_pos_dep': self.case_depth,
            'Ctrl_pos_dep': self.control_depth,
            'Target_gene': self.target_gene,
            'Scale_ctrl': self.control_score,
            'rsID': self.rs,
            'TMB_Type': self.tmb,
            'Bl_MutType': self.bl_muttype,
            'filter_0.01': self.filter_001,
            'filter_0.05': self.filter_005,
            'Interpret': self.interpret,
            '是否是库内点': self.in_tianhua,
            '入库时间': self.tianhua_time,
            '同一坐标变异': self.same_position_variant,
            '龙舌兰文献': self.tequila_literature,
            '龙舌兰分级': self.tequila_class,
            '龙舌兰ACMG等级': self.tequila_acmg,
            'Clinvar': self.clinvar,
            '引物设计': self.primer_design
        }

        info.update(tags)

        return info



if __name__ == '__main__':
    head = HEAD()
    head.gene = 'EGFR'
    head_info = head.update_head()
    with open('head.xls', 'w') as fw:
        fw.write('\t'.join(head_info.keys()) + '\n')
        for key in head_info:
            fw.write(head_info[key] + '\t')
        
        