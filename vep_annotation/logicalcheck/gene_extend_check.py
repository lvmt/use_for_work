#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   gene_extend_check.py
@Time    :   2021/12/13 15:44:02
@Author  :   lvmt 
@Version :   1.0
'''

# here put the import lib
'''
变异发生后，氨基酸序列是否延长
主要争对区域：frameshift
'''



from os import pardir
import re  


class GeneExtendCheck:
    
    def __init__(self, args, tran, func, phgvs):
        self.gene_extend_database = args['gene_extend_confg']
        self.tran = tran
        self.func = func  
        self.phgvs = phgvs 


    def get_gene_extend_info(self):
        gene_extend_info = {}
        with open(self.gene_extend_database, 'r') as fr:
            for line in fr:
                if line.startswith('#'):
                    continue
                linelist = line.strip('\n').split('\t')
                tran = linelist[1]
                tran_length = linelist[2]
                gene_extend_info[tran] = tran_length
            
        return gene_extend_info


    def get_newlength(self):
        '''判断发生fs变异后，新蛋白质的长度
        '''
        # pattern = re.compile(r'p.[A-Z]+(\d+)[A-Za-z\*]+(\d+)') # p.Q118Kfs*4 
        pattern = re.compile(r'p.[A-Z]+(\d+)[A-Za-z\*]+(\d*)')  #  特例：p.Y370*
        if re.search(pattern, self.phgvs):
            # pos1, pos2 = re.search(pattern, self.phgvs).groups()  # 特例：p.Y370*
            pos_list = re.search(pattern, self.phgvs).groups()  # 特例：p.Y370*
            return sum(map(int, filter(str, pos_list))) - 1 
        return 'compile wrong'
 

    def start(self):
        gene_extend_info = self.get_gene_extend_info()
        if self.func == 'frameshift':
            origin_length = gene_extend_info.get(self.tran, 20000)  # 防止库文件没有该信息
            new_length = self.get_newlength()
            if new_length > int(origin_length):
                return 'Y'
            else:
                return 'N'
        return 'Nan'



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='just a test')
    parser.add_argument('--infile')
    parser.add_argument('--gene_extend_confg')

    args = vars(parser.parse_args())

    with open(args['infile'], 'r') as fr, open('test.add_gene_extend', 'w') as fw:
        for line in fr:
            if line.startswith('Scale'):
                fw.write(line)
                continue
            linelist = line.strip('\n').split('\t')
            tran = linelist[37]
            func = linelist[9]
            phgvs = linelist[3] 

            linelist[12] = GeneExtendCheck(args, tran, func, phgvs).start()

            fw.write('{}\n'.format('\t'.join(linelist)))
             
