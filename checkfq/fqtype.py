#!/usr/bin/env python3
#-*- coding:utf-8 -*-

"""
0. 根据sample_list, 得到包含下机路径的列表
1. 根据fq文件首行，进行判断，返回nova 还是x-ten
2. 如果测序是nova，输出不是nova的
3. 如果测序是x-ten，输出不是x-ten的
4. need python3
"""

import argparse

class CheckType(object):

    def __init__(self, args):
        self.infile = args['infile']
        self.type = args['type']
        self.info_list = self.get_info_list()

    def get_info_list(self):
        info_list = []
        with open(self.infile) as f:
            for line in f:
                if line.startswith('#'):
                    # title = line.strip().split('\t')
                    pass
                else:
                    linelist = line.strip().split('\t')
                    fq_path = "/".join((linelist[6], linelist[3], linelist[3])) + "_L" + linelist[0] + "_1.fq.gz"
                    info_list.append(fq_path)
        
        return info_list
 
    def read_title(self, fq):
        try:
            if fq.endswith('gz'):
                import gzip
                return gzip.open(fq).readline().strip()
            else:
                return open(fq).readline().strip()
        except Exception as e:
            print(e)
            exit("func wrong: infile is wrong")

    def print_info(self, stat_info):
        """
        stat info list
        """
        if stat_info:
            cell_set = set(stat_info)
            n = 1 
            print("\033[31m 下面路径格式存在问题\033[0m")
            for item in cell_set:
                print(n, item, end="\n")
                n += 1
        else:
            print('\033[31m 数据格式没有问题\033[0m')

    def check(self):
        nova_info = []
        x_ten_info = []

        for fq in self.info_list:
            title = self.read_title(fq)
            fq = "/".join((fq.split('/', 6))[0:-1])
            if title.startswith(b"@ST"):     # x-ten format 
                x_ten_info.append(':'.join(('x-ten', fq)))
            else:
                nova_info.append(':'.join(('nova', fq)))
    
        if self.type.lower() == "nova":
           self.print_info(x_ten_info)
        elif self.type.lower() == "x-ten":
            self.print_info(nova_info)
        else:
            print("type: x-ten or nova")


def main():

    check = CheckType(args)
    check.check()


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="check fq type, nova or x-ten")
    parser.add_argument('infile', help="sample_list")
    parser.add_argument('type', help="fq, ", choices=['nova', 'x-ten'])
    args = vars(parser.parse_args())

    main()




