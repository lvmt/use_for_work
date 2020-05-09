#!/usr/bin/env python
#-*- coding:utf-8 -*-

"""
文库数据量统计
"""

import os
import subprocess

class StatG(object):

    def __init__(self, args):

        self.infile = args['infile']
        self.statfile = args['statfile']

    def stat(self):
        with open(self.infile, 'r') as fr, open(self.statfile, 'w') as fw:
            for line in fr:
                line = line.strip()
                if line.startswith('#'):
                    fw.write('{line}\t数据量大小\n'.format(**locals()))
                    continue
                linelist = line.strip().split('\t')
                lib = linelist[58]
                index = linelist[44]
                path = linelist[66] 
                path = self.get_path(lib, index, path)

                if path != "None":
                    path_size = self.get_path_size(path)
                else:
                    path_size = "None"
                fw.write("{line}\t{path_size}\n".format(**locals()))
                # print(path, path_size)

    def get_path_size(self, path):
        stat_file = str(path) + "/*2.adap.stat"
        stat_file = subprocess.getstatusoutput("ls {stat_file}".format(**locals()))[1]
        with open(stat_file, 'r') as fr:
            for line in fr:
                if line.startswith("total_reads"):
                    path_size = line.strip().split(":")[-1].strip()
                    path_size = int(path_size) * 300  / 1000000000 
                else:
                    pass
        return path_size

    def get_path(self, lib, index, path):

        lib = lib.strip().replace(";", "-")
        index = index.strip().replace(";", "-")
        path = path.strip()
        
        newpath1 = "/".join((path, lib))
        newpath2 = path + "/" + lib + "-" + index

        if os.path.exists(newpath1):
            return newpath1
        elif os.path.exists(newpath2):
            return newpath2
        else:
            # print("该路径有问题:\n{path} {lib} {index}".format(**locals()))
            return "None"
            


def main():

    demo = StatG(args)
    demo.stat()


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--infile", help="输入文件")
    parser.add_argument("--statfile", help="输出统计文件")

    args = vars(parser.parse_args())

    main()



