#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   splice_affect_check.py
@Time    :   2021/12/13 15:45:13
@Author  :   lvmt 
@Version :   1.0
'''

# here put the import lib
'''
变异是否影响剪切
check标签
 - splice-3
 - splice-5
 - span

根据梳理出来的逻辑直接转换成代码
其他标签输出Nan
'''


from os import pardir
import re
import yaml



#############################################
############## vep中各个变异示例 #############
# del: c.492_494del
# ins: c.2310_2311insT / c.583_584insCCCCGC
# delins: c.35_36delinsCT
# snp: c.126T>C

#################################################
############### 可能存在的情况枚举 ################

## splice: del 类
# pattern1 = re.compile(r'c.(\d+)\+1_(\d+)\+\d+del')      # c.111+1_111+5del
# pattern2 = re.compile(r'c.(\d+)\+2_(\d+)\+\d+del')      # c.111+2_111+5del
# pattern3 = re.compile(r'c.(\d+)-\d+_(\d+)-1del')        # c.222-5_222-1del
# pattern4 = re.compile(r'c.(\d+)-\d+_(\d+)-2del')        # c.222-5_222-2del

pattern1 = re.compile(r'c.(\d+)\+1[_]*(\d*)\+*\d*del')      #c.111+1_111+5del, c.111+1del
pattern2 = re.compile(r'c.(\d+)\+2[_]*(\d*)\+*\d*del')      # c.111+2_111+5del, c.111+2del 
pattern3 = re.compile(r'c.(\d*)-*\d*[_]*(\d+)-1del')        # c.222-5_222-1del, c.222-1del
pattern4 = re.compile(r'c.(\d*)-*\d*[_]*(\d+)-2del')        # c.222-5_222-2del, c.222-2del
# pattern13 = re.compile(r'c.(\d+)\+1del')                  # c.111+1del
# pattern14 = re.compile(r'c.(\d+)\+2del')                  # c.111+2del
# pattern15 = re.compile(r'c.(\d+)-1del')                   # c.222-1del
# pattern16 = re.compile(r'c.(\d+)-2del')                   # c.222-2del


## splice: delins 类
# 只要是发生在splice区域的delins，都会影响剪切
# pattern2 = re.compile(r'c.(\d+)\+1_(\d+)\+(\d)delins')   # c.111+1_111+5delins
# pattern5 = re.compile(r'c.(\d+)-(\d+)_(\d)-1delins')     # c.222-5_222-1delins
# pattern17 = re.compile(r'c.(\d+)\+1delins')              # c.111+1delins
# pattern20 = re.compile(r'c.(\d+)-1delins')               # c.222-1delins
# pattern21 = re.compile(r'c.(\d+)\+2delins')              # c.111+2delins
# pattern21 = re.compile(r'c.(\d+)-2delins')               # c.222-2delins


## splice: ins类
pattern5 = re.compile(r'c.(\d+)_(\d+)\+1ins')     # c.111_111+1ins
pattern6 = re.compile(r'c.(\d+)\+1_(\d+)\+2ins')  # c.111+1_111+2ins
pattern7 = re.compile(r'c.(\d+)-2_(\d+)-1ins')    # c.222-2_222-1ins
pattern8 = re.compile(r'c.(\d+)-1_(\d+)ins')      # c.222-1_222ins 
pattern9 = re.compile(r'c.(\d+)\+1dup')           # c.222+1dup
pattern10 = re.compile(r'c.(\d+)\+2dup')          # c.222+2dup
pattern11 = re.compile(r'c.(\d+)-1dup')           # c.333-1dup
pattern12 = re.compile(r'c.(\d+)-2dup')           # c.333-2dup


## splice: snp 类
pattern13 = re.compile(r'c.(\d+)\+2[A-Z]>')  # c.111+2T>C
pattern14 = re.compile(r'c.(\d+)\+1[A-Z]>')  # c.111+1G>A
pattern15 = re.compile(r'c.(\d+)-2[A-Z]>')   # c.111-2A>C
pattern16 = re.compile(r'c.(\d+)-1[A-Z]>')   # c.111-1G>A


#################################################################
################### span类型 可能存在情况枚举 #####################

## span: delins 
pattern17 = re.compile(r'c.(\d+)-\d+_(\d+)delins')    # c.111-4_115delins (span)
pattern18 = re.compile(r'c.(\d+)_(\d+)\+\d+delins')   # c.88_95+3delins   (span)
 

## span: del （pattern19和pattern22存在歧义：无法确定右端是位于紧邻exon,还是下下个exon）
pattern19 = re.compile(r'c.(\d+)-\d+_(\d+)del')           # c.111-4_115del (span)
pattern20 = re.compile(r'c.(\d+)-\d+_(\d+)\+\d+del')      # c.111-4_200+5del  (span,跨越exon)
pattern21 = re.compile(r'c.(\d+)_(\d+)\+\d+del')          # c.98_110+4del 
pattern22 = re.compile(r'c.(\d+)\+\d+_(\d+)del')          # c.111+30_120del
pattern23 = re.compile(r'c.(\d+)\+\d+_(\d+)\+\d+del')     # c.110+3_300+5del (span, 跨越exon)


splice_info = {
    'del': {
        'pattern1': pattern1,
        'pattern2': pattern2,
        'pattern3': pattern3,
        'pattern4': pattern4
    },

    'ins':  {
        'pattern5': pattern5,
        'pattern6': pattern6,
        'pattern7': pattern7,
        'pattern8': pattern8,
        'pattern9': pattern9,
        'pattern10': pattern10,
        'pattern11': pattern11,
        'pattern12': pattern12,
    },

    'snp': {
        'pattern13': pattern13,
        'pattern14': pattern14,
        'pattern15': pattern15,
        'pattern16': pattern16
    }
}


span_info = {
    'delins': {
        'pattern17': pattern17,
        'pattern18': pattern18,
    },

    'del': {
        'pattern19': pattern19,
        'pattern20': pattern20,
        'pattern21': pattern21,
        'pattern22': pattern22,
        'pattern23': pattern23
    }
}


affect_pattern_list = ['pattern9', 'pattern11', 'pattern13', 'pattern14', 'pattern15', 'pattern16', 'pattern20', 'pattern23']
not_affect_pattern_list = ['pattern10', 'pattern12']



class SpliceAffectCheck:

    def __init__(self, args, func, chgvs, flank, strand):
        self.splice_pattern_yaml = args['splice_affect_config']
        self.func = func 
        self.chgvs = chgvs 
        self.flank = flank
        self.strand = strand 


    def not_UTR_variant(self, chgvs):
        '''判断该变异是否是UTR区域
        对应UTR区域的splice和span,对于剪切没有影响
        '''
        if re.search(r'c\.-\d', chgvs) or re.search(r'c\.\*\d', chgvs) \
            or re.search(r'_-(\d+)', self.chgvs) or re.search(r'_\*(\d_)', self.chgvs):
            return False
        return True


    def get_target_pattern(self, pattern_info):
        '''
        pattern_info
        {
            'pattern1': pattern1, 'pattern2': pattern2
        }
        '''
        for pattern in pattern_info:
            if re.search(pattern_info[pattern], self.chgvs):
                return pattern

        return 'None'  # 暂时不知道是否存在额外情况


    def get_muttype(self):
        if 'delins' in self.chgvs:
            return 'delins'
        elif 'del' in self.chgvs:
            return 'del'
        elif 'ins' in self.chgvs or 'dup' in self.chgvs:
            return 'ins'
        else:
            return 'snp'


    def get_pattern_config(self):
        pattern_config = yaml.load(open(self.splice_pattern_yaml))
        return pattern_config


    def handle_splice_func(self, muttype, pattern_config):
        '''muttype根据chgvs内容判断
        '''
        if 'delins' in muttype:  # splice中的delins，全都影响剪切
            return 'Y'

        pattern_info = splice_info[muttype] 
        target_pattern = self.get_target_pattern(pattern_info)
        # print('测试:', target_pattern)

        ## 
        if target_pattern in affect_pattern_list:
            return 'Y'
        elif target_pattern in not_affect_pattern_list:
            return 'N'
            
        left, right = self.flank.upper().split('.')
        ins = ''
        if muttype == 'ins':
            ins = re.search(r'ins([a-zA-Z]+)', self.chgvs).group(1)

        # 改变异位点本身的侧翼序列及其插入序列
        base_info = {
            'left': left,
            'right': right,
            'ins': ins
        }

        key = list(pattern_config[target_pattern][self.strand].keys())[0]  # left, right, ins
        value = list(pattern_config[target_pattern][self.strand].values())[0] # 序列信息 ^T, $A等

        if re.search(r'{value}'.format(**locals()), base_info[key]):
            return 'N'
        else:
            return 'Y'


    def handle_span_func(self, muttype, pattern_config):

        pattern_info = span_info[muttype]
        target_pattern = self.get_target_pattern(pattern_info)  # 字符串
        # print('测试:', target_pattern)

        if target_pattern in affect_pattern_list:
            return 'Y'
        elif target_pattern in not_affect_pattern_list:
            return 'N'

        left, right = self.flank.upper().split('.')
        ins = ''
        if muttype == 'delins':
            ins = re.search(r'delins([a-zA-Z]+)', self.chgvs).group(1)

        base_info = {
            'left': left,
            'right': right,
            'ins': ins
        }

        ## 优先判断是否为3的倍数 
        pos1, pos2 = re.search(pattern_info[target_pattern], self.chgvs).groups()
        if not (abs(int(pos1) - int(pos2)) + 1) % 3 == 0:
            return 'Y'
        else:
            key = list(pattern_config[target_pattern][self.strand].keys())[0]  # left, right, ins
            value = list(pattern_config[target_pattern][self.strand].values())[0] # 序列信息 ^T, $A等

            if re.search(r'{value}'.format(**locals()), base_info[key]):
                return 'N'
            else:
                return 'Y'

    
    def start(self):
        '''
        优先判断该span类型位于5'UTR或者3'UTR
        c.-20_16delinsCGGCGGCGGAGGAGGCGGCGACCGAGAAGATGCCCGCCCTGCG
        '''
        muttype = self.get_muttype()
        pattern_config = self.get_pattern_config()

        if self.func.startswith('splice') and self.not_UTR_variant(self.chgvs):
            return self.handle_splice_func(muttype, pattern_config)
        elif self.func.startswith('span') and self.not_UTR_variant(self.chgvs):
            return self.handle_span_func(muttype, pattern_config)
        else:
            return 'Nan'



if __name__ == '__main__':

    # func = 'splice-3'
    # chgvs = 'c.111+1_111+5del'
    # flank = 'Ac.ac'
    # strand = '-'
    splice_config = 'splice.yaml'

    ## splice delins 测试
    # print(SpliceAffectCheck('splice-3', 'c.222-8_222-2delins', 'aa.tt', '-', splice_config).start())  # affect

    # pattern1 测试
    # print(SpliceAffectCheck(splice_config, 'splice-3', 'c.111+1del', 'Ac.tc', '-').start())  # 没
    
    # print(SpliceAffectCheck('splice-3', 'c.111+1_111+5del', 'gt.ac', '-', splice_config).start())  # 有
    # print(SpliceAffectCheck('splice-3', 'c.111+1_111+5del', 'gt.ac', '+', splice_config).start())  # 有
    # print(SpliceAffectCheck('splice-3', 'c.111+1_111+5del', 'gt.gt', '+', splice_config).start())  # 没

    ## pattern2 测试
    # print(SpliceAffectCheck(splice_config, 'splice-3', 'c.111+2del', 'Ac.tc', '-').start())  # 没
    # print(SpliceAffectCheck('splice-3', 'c.111+2_111+5del', 'gt.gt', '+', splice_config).start())  # 有
    # print(SpliceAffectCheck('splice-3', 'c.111+2_111+5del', 'gt.tt', '+', splice_config).start())  # 没
    # print(SpliceAffectCheck('splice-3', 'c.111+2_111+5del', 'at.tt', '-', splice_config).start())  # 有
    # print(SpliceAffectCheck('splice-3', 'c.111+2_111+5del', 'aa.tt', '-', splice_config).start())  # 没
    
    ## pattern3 测试
    # print(SpliceAffectCheck(splice_config, 'splice-3', 'c.222-1del', 'Ac.tc', '-').start())  # 没
    # print(SpliceAffectCheck('splice-3', 'c.222-5_222-1del', 'aa.ag', '+', splice_config).start())  # 有
    # print(SpliceAffectCheck('splice-3', 'c.222-5_222-1del', 'ag.ag', '+', splice_config).start())  # 没
    # print(SpliceAffectCheck('splice-3', 'c.222-5_222-1del', 'ag.ag', '-', splice_config).start())  # 有
    # print(SpliceAffectCheck('splice-3', 'c.222-5_222-1del', 'ag.ct', '-', splice_config).start())  # 没

    ## pattern4 测试
    print(SpliceAffectCheck(splice_config, 'splice-3', 'c.222-2del', 'Ac.tc', '-').start())  # 没
    # print(SpliceAffectCheck('splice-3', 'c.222-8_222-2del', 'ag.ct', '+', splice_config).start())  # 有
    # print(SpliceAffectCheck('splice-3', 'c.222-8_222-2del', 'aa.ct', '+', splice_config).start())  # 没
    # print(SpliceAffectCheck('splice-3', 'c.222-8_222-2del', 'aa.ct', '-', splice_config).start())  # 有
    # print(SpliceAffectCheck('splice-3', 'c.222-8_222-2del', 'aa.tt', '-', splice_config).start())  # 没


    ## pattern5 测试
    # print(SpliceAffectCheck('splice-3', 'c.111_111+1insTT', 'aa.tt', '+', splice_config).start())  # 有
    # print(SpliceAffectCheck('splice-3', 'c.111_111+1insGT', 'aa.tt', '+', splice_config).start())  # 没
    # print(SpliceAffectCheck('splice-3', 'c.111_111+1insGTA', 'aa.tt', '+', splice_config).start())  # 没
    # print(SpliceAffectCheck('splice-3', 'c.111_111+1insTA', 'aa.tt', '-', splice_config).start())  # 有
    # print(SpliceAffectCheck('splice-3', 'c.111_111+1insGT', 'aa.tt', '-', splice_config).start())  # 没


    ## pattern6 测试
    # print(SpliceAffectCheck('splice-3', 'c.111+1_111+2insAA', 'aa.tt', '-', splice_config).start())  # 有
    # print(SpliceAffectCheck('splice-3', 'c.111+1_111+2insTT', 'aa.tt', '-', splice_config).start())  # 没

    ## pattern7 测试
    # print(SpliceAffectCheck('splice-3', ' c.222-2_222-1insAAT', 'aa.tt', '-', splice_config).start())  # 有
    # print(SpliceAffectCheck('splice-3', ' c.222-2_222-1insAA', 'aa.tt', '-', splice_config).start())  # 没

    ## pattern8 测试
    # print(SpliceAffectCheck('splice-3', 'c.222-1_222insAA', 'aa.tt', '-', splice_config).start())  # 有
    # print(SpliceAffectCheck('splice-3', 'c.222-1_222insAG', 'aa.tt', '-', splice_config).start())  # 没

    ## pattern9 测试
    # print(SpliceAffectCheck('splice-3', 'c.222+1dup', 'aa.tt', '-', splice_config).start())  # 有

    ## pattern10 测试
    # print(SpliceAffectCheck('splice-3', 'c.222+2dup', 'aa.tt', '-', splice_config).start())  # 肯定不影响

    ## pattern11 测试
    # print(SpliceAffectCheck('splice-3', 'c.333-1dup', 'aa.tt', '-', splice_config).start())  # 肯定影响
    # print(SpliceAffectCheck('splice-3', 'c.333-1dup', 'aa.tt', '+', splice_config).start())  # 肯定影响

    ## pattern12 测试
    # print(SpliceAffectCheck('splice-3', 'c.333-2dup', 'aa.tt', '+', splice_config).start())  # 肯定不影响

    ## pattern13 ~ 16 测试
    # print(SpliceAffectCheck('splice-3', 'c.111+2T>C', 'aa.tt', '+', splice_config).start())  # 肯定影响
    # print(SpliceAffectCheck('splice-3', 'c.111+1T>C', 'aa.tt', '+', splice_config).start())  # 肯定影响
    # print(SpliceAffectCheck('splice-3', 'c.111-2A>C', 'aa.tt', '+', splice_config).start())  # 肯定影响
    # print(SpliceAffectCheck('splice-3', 'c.111-1G>C', 'aa.tt', '+', splice_config).start())  # 肯定影响


    ## pattern17 测试 
    # print(SpliceAffectCheck('span', 'c.111-4_115delinsATCG', 'aa.tt', '+', splice_config).start())  # 不等于3，影响剪切
    # print(SpliceAffectCheck('span', 'c.111-8_115delinsATCG', 'aa.tt', '+', splice_config).start())  # 不等于3，影响剪切
    # print(SpliceAffectCheck('span', 'c.111-8_116delinsATCG', 'aa.tt', '+', splice_config).start())    # 存在影响
    # print(SpliceAffectCheck('span', 'c.111-8_116delinsATCAG', 'aa.tt', '+', splice_config).start())   # 没有影响
    # print(SpliceAffectCheck('span', 'c.111-8_116delinsATCAG', 'aa.tt', '-', splice_config).start())   # 没有影响


    ## pattern18 测试 
    # print(SpliceAffectCheck('span', 'c.88_95+3delinsATCAG', 'aa.tt', '-', splice_config).start())   # 不等于3，影响剪切
    # print(SpliceAffectCheck('span', 'c.88_96+3delinsATCAG', 'aa.tt', '-', splice_config).start())   # 存在影响
    # print(SpliceAffectCheck('span', 'c.88_96+3delinsGTCAG', 'aa.tt', '-', splice_config).start())   # 没有影响
    # print(SpliceAffectCheck('span', 'c.88_96+3delinsGTCAG', 'aa.tt', '+', splice_config).start())   # 没有影响


    ## pattern19 测试 
    # print(SpliceAffectCheck('span', 'c.111-4_115del', 'aa.tt', '+', splice_config).start())   # 不等于3，影响剪切
    # print(SpliceAffectCheck('span', 'c.111-4_116del', 'aa.tt', '+', splice_config).start())   # 存在影响
    # print(SpliceAffectCheck('span', 'c.111-4_116del', 'ag.tt', '+', splice_config).start())   # 没有影响
    # print(SpliceAffectCheck('span', 'c.111-4_116del', 'ag.tt', '-', splice_config).start())   # 存在影响
    # print(SpliceAffectCheck('span', 'c.111-4_116del', 'ag.ct', '-', splice_config).start())   # 没有影响

    ## pattern20 测试
    # print(SpliceAffectCheck('span', 'c.111-4_200+5del', 'ag.ct', '-', splice_config).start())   # 肯定影响

    ## pattern21 测试 
    # print(SpliceAffectCheck('span', 'c.98_110+4del', 'ag.ct', '-', splice_config).start())   # 不等于3，影响剪切
    # print(SpliceAffectCheck('span', 'c.98_112+3del', 'ag.ct', '+', splice_config).start())   # 存在影响
    # print(SpliceAffectCheck('span', 'c.98_112+3del', 'ag.gt', '+', splice_config).start())   # 没有影响
    # print(SpliceAffectCheck('span', 'c.98_112+3del', 'ag.ac', '-', splice_config).start())   # 存在影响
    # print(SpliceAffectCheck('span', 'c.98_112+3del', 'ac.ac', '-', splice_config).start())   # 没有影响

    ## pattern22 测试 
    # print(SpliceAffectCheck('span', 'c.111+30_120del', 'ac.ac', '-', splice_config).start())   # 不等于3，影响剪切
    # print(SpliceAffectCheck('span', 'c.111+30_119del', 'ag.ac', '+', splice_config).start())   # 没有影响
    # print(SpliceAffectCheck('span', 'c.111+30_119del', 'ag.ac', '-', splice_config).start())   # 存在影响
    # print(SpliceAffectCheck('span', 'c.111+30_119del', 'ag.ct', '-', splice_config).start())   # 没有影响

    ## pattern23 测试 
    # print(SpliceAffectCheck('span', ' c.110+3_300+5del', 'ag.ct', '-', splice_config).start())   # 肯定影响

    ## 其实func测试
    # print(SpliceAffectCheck('missense', ' c.110A>T', 'ag.ct', '-', splice_config).start())   # no affect

