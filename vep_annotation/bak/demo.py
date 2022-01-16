#!/usr/bin/python

import re 
import argparse
import utils

class ROW:
    def __init__(self):
        self.name = '*'
        self.age = '*'
        self.sex = '*'

    def update_head(self, **args):
        self.__dict__.update(args)      
        return self


tmp = ROW()
tmp.name = 'mike'
tmp.age = 23

tag = {
    'class': 'A',
    'like': 'apple'
}

info = tmp.update_head(**tag)

print(info)
print(info.__dict__)
    

### test oneletter acid
# hgvsp_list = ['p.Thr26Met', 'p.Thr26delinsIleSer', 'p.Thr26ArgfsTer167', 
#               'p.Thr26ValfsTer165', 'p.Thr26del', 'p.Thr26Gly', 'p.Thr26ArgfsTer166']

# for hgvsp in hgvsp_list:
#     print(utils.get_oneletter_hgvsp(hgvsp))


### test reverse_com
# base = 'ATG'
# print(utils.reverse_complement(base))


### test snv vcf info 
# tmp = 'BASE_DLY=0;DP4=1448,90,26,3;CtrlDP4=463,29,5,0;DUPLEX=9,4,8,0'

# case_pattern = re.compile(r'DP4=(\d+),(\d+),(\d+),(\d);?')
# control_pattern = re.compile(r'CtrlDP4=(\d+),(\d+),(\d+),(\d);?')

# case = re.search(case_pattern, tmp)
# control = re.search(control_pattern, tmp)

# ll = []

# if case:
#     case = case.group(1)
#     ll.append(case)
# if control:
#     control = control.groups()
#     ll.append(control)

# print(ll)

### test modify_pos_ref_alt
# print('ins', utils.modify_pos_ref_alt(1295028, 'G', 'GAAA'))
# print('del', utils.modify_pos_ref_alt(1295028, 'GTGG', 'G'))
# print('delins', utils.modify_pos_ref_alt(1295028, 'GT', 'CC'))
# print('delins', utils.modify_pos_ref_alt(1295028, 'G', 'CC'))

