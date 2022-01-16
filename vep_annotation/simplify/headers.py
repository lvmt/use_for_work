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
        self.phgvs_shoft = '*'
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
        self.tert = '*'
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
            'Hugo_Symbol': self.gene,
            'HGVSc': self.chgvs,
            'HGVSp_Short': self.phgvs_shoft,
            'HGVSp': self.phgvs,
            'EXON': self.exon_id,
            'Case_var_freq': self.case_var_freq,
            'VEP_Function': self.vep_function,
            'Consequence': self.vep_simple_function,
            'BGI_Function': self.vep2bgicg_function,
            'BGI_NewFunction': self.newfunction,
            'BGI_End_Exon_Check': self.end_exon_check,
            'BGI_Gene_Extend_Check': self.gene_extend_check,
            'BGI_Splice_Affect_Check': self.splice_affect_check,
            'BGI_Uniport_Position(s)': self.uniport_position,
            'BGI_Uniport_feature_key': self.uniport_featurn,
            'BGI_Uniport_description': self.uniport_description,
            'Cosmic': self.cosmic,
            # 'Tert': self.tert,
            # 'Cosmic2': self.cosmic2,
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
            'Chromosome': self.chr,
            'Start_Position': self.start,
            'End_Position': self.end,
            'Reference_Allele': self.ref,
            'Allele': self.alt,
            'MapLoc': self.maploc,
            'BGI_Variant_Type': self.muttype,
            'Gtype': self.Gtype,
            'Codons': self.genotype,
            'Transcript_ID': self.transcript,
            'ENSP': self.protein,
            'Strand': self.strand,
            'FuncRegion': self.funcregion,
            'Flank': self.flank,
            'Case_pos_dep': self.case_depth,
            'Ctrl_pos_dep': self.control_depth,
            'Target_gene': self.target_gene,
            'Scale_ctrl': self.control_score,
            'rsID': self.rs,
            'TMB_Type': self.tmb,
            'Variant_Type': self.bl_muttype,
            'Local_AF': '*',
            'Local_Hom_AF': '*',
            'GMAF': self.af,
            'EAS_AF': self.eas_af,
            'AMR_AF': self.amr_af,
            'AFR_AF': self.afr_af,
            'EUR_AF': self.eur_af,
            'SAS_AF': self.sas_af,
            'ESP6500': self.ea_af,
            'gnomAD_AF': self.gnomad_af,
            'gnomAD_EAS_AF': self.gnomad_eas_af,
            'gnomAD_EAS_Hom_AF': self.gnomad_afr_af,
            'gnomAD_AFR_AF': self.gnomad_afr_af,
            'gnomAD_AFR_Hom_AF': self.gnomad_afr_af,
            'gnomAD_AMR_AF': self.gnomad_amr_af,
            'gnomAD_AMR_Hom_AF': '*',
            'gnomAD_FIN_AF': self.gnomad_fin_af,
            'gnomAD_FIN_Hom_AF': '*',
            'gnomAD_NFE_AF': self.gnomad_nfe_af,
            'gnomAD_NFE_Hom_AF': '*',
            'gnomAD_SAS_AF': self.gnomad_sas_af,
            'gnomAD_SAS_Hom_AF': '*',
            'filter_0.01': self.filter_001,
            'filter_0.05': self.filter_005,
            'Interpret': self.interpret,
            'InDatabase_check': self.in_tianhua,
            'InDatabase_time': self.tianhua_time,
            'Same_coordinate_variation ': self.same_position_variant,
            'tequila_literature': self.tequila_literature,
            'tequila_class': self.tequila_class,
            'tequila_acmg': self.tequila_acmg,
            'Clinvar': self.clinvar,
            'Primer_design': self.primer_design
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
        
        