#!/usr/bin/env python
# -*- coding:utf-8 -*-


"""
对比广妇儿项目修改流程结果的正确性
"""


class Diff(object):

    def __init__(self, args):
        self.file1 = args['file1']
        self.file2 = args['file2']

    @staticmethod
    def get_list(filename):
        str_list = []
        with open(filename, 'r') as f:
            for line in f:
                if line.startswith('Priority'):
                    continue
                linelist = line.strip().split('\t')
                new_str = "|".join(linelist[0:6])
                str_list.append(new_str)
        return str_list

    def diff(self):
        diff_result = []
        file1_str_list = Diff.get_list(self.file1)
        file2_str_list = Diff.get_list(self.file2)

        for item in file1_str_list:
            if item in file2_str_list:
                pass
            else:
                diff_result.append((item, "not_file2"))

        for item in file2_str_list:
            if item in file1_str_list:
                pass
            else:
                diff_result.append(('not_file1', item))

        print(len(diff_result))
        print(diff_result)


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('file1')
    parser.add_argument('file2')
    args = vars(parser.parse_args())

    diff = Diff(args)
    diff.diff()


if __name__ == "__main__":
    main()
