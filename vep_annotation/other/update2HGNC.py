#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   update2HGNC.py
@Time    :   2022/01/05 15:58:12
@Author  :   lvmt 
@Version :   1.0
'''

# here put the import lib
# 更正基因名为标准HGNC名称


class ModifyName:

    def __init__(self, args):
        self.infile = args['infile']
        self.result = args['result']
        self.hgnc = args['hgnc_config']

    
    def get_head_index(self, head):
        head_index = {}
        if not isinstance(head, list):
            head = head.strip('\n').split('\t')

        for index,item in enumerate(head):
            head_index[item.lower()] = index 

        return head_index


    def get_name_relation(self):
        relation_info = {}
        with open(self.hgnc, 'r') as fr:
            for line in fr:
                linelist = line.strip('\n').split('\t')
                hgnc_name = linelist[1]
                ncbi_name = linelist[2]
                relation_info[ncbi_name] = hgnc_name

        return relation_info


    def start(self):
        relation_info = self.get_name_relation()

        with open(self.infile, 'r') as fr, open(self.result, 'w') as fw:
            for line in fr:
                if line.startswith('Scale'):
                    fw.write(line)
                    head_index = self.get_head_index(line)
                    continue
                linelist = line.strip('\n').split('\t')
                gene = linelist[head_index['hugo_symbol']]
                linelist[head_index['hugo_symbol']] = relation_info.get(gene, gene)  # 如何不在,还是使用自己
                fw.write('{}\n'.format('\t'.join(linelist)))



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='修改基因名为HGNC标准名称')
    parser.add_argument('--infile', help='输入')
    parser.add_argument('--result', help='输出')
    parser.add_argument('--hgnc_config', help='基因名对应关系配置文件')

    args = vars(parser.parse_args())

    ModifyName(args).start()
    

