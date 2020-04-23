#!/usr/bin/env python
#-*- coding:utf-8 -*-

"""
随机提起1000个位点
"""

import random
import argparse

class Extract(object):

    def __init__(self, args):
        self.infile = args['infile']
        self.out = args['out']
    
    @staticmethod
    def safe_open(infile):
        try:
            if infile.endswith('gz'):
                import gzip
                return gzip.open(infile)
            else:
                return open(infile)
        except Exception as e:
            print(e)
            exit("something wrong") 

    def write_file(self):
        f = Extract.safe_open(self.infile)
        o = open(self.out, 'w')

        ## write title
        for line in f:
            if line.startswith('#'):
                o.write(line)
            else:
                break
            
        file_list = f.readlines()
        random_set = sorted(random.sample(range(2, 18000), 1000))  
        print(len(random_set))
        for index, line in enumerate(file_list):
            if index in random_set:
                o.write(line)
        
def main():
    e = Extract(args)
    e.write_file()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile')
    parser.add_argument('out')
    args = vars(parser.parse_args())

    main()