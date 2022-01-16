#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   cismut_filter_normal.py
@Time    :   2022/01/06 15:41:32
@Author  :   lvmt 
@Version :   1.0
'''

# here put the import lib
# 顺势结果,如果control的支持数>=10，且变异类型不是snv,则舍弃

import re  



def get_head_index(head):
    head_info = {}
    if not isinstance(head, list):
        head = head.strip('\n').split('\t')
    
    for index,item in enumerate(head):
        head_info[item.lower()] = index

    return head_info


def start(infile, outfile):
    with open(infile, 'r') as fr, open(outfile, 'w') as fw:
        for line in fr:
            if line.startswith('Scale'):
                fw.write(line)
                head_info = get_head_index(line)
                continue
            linelist = line.strip('\n').split('\t')
            muttype = linelist[head_info['bgi_variant_type']]
            normal_var_read = linelist[head_info['ctrl_var_readsnum']]

            if re.search(r'(\d+)', normal_var_read):
                normal_alt_read = re.search(r'(\d+)', normal_var_read).group(1)
                if float(normal_alt_read) >= 10 and not muttype == 'snv':
                    continue
            linelist[0] = 'Yes'
            fw.write('{}\n'.format('\t'.join(linelist)))



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='过滤对照样本变异支持数')
    parser.add_argument('--infile')
    parser.add_argument('--result')

    args = vars(parser.parse_args())

    start(args['infile'], args['result'])
    
