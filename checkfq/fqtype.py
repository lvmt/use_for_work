#!/usr/bin/env python
#-*- coding:utf-8 -*-

class CheckDataType(object):

    def __init__(self, args):
        self.infile = args['infile']   # sample_list

    @staticmethod
    def get_title(infile):
        try:
            if not infile.endswith('gz'):
                return open(infile).readline()
            else:
                import gzip
                return gzip.open(infile).readline()
        except Exception as e:
            print(e)
        
    def get_sample_dict(self):
        with open(self.infile) as f:

        
        
# wgs : 
#     snp : 395670
#     indel : 76293


    

