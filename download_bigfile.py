#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   download_bigfile.py
@Time    :   2021/08/27 17:50:25
@Author  :   lvmt
@Version :   1.0
'''

'''
对大文件进行切割下载
文件总大小，一次下载大小，分割文件前缀，下载文件路径
'''

# here put the import lib


class Download:

    def __init__(self, args):
        self.block_size = args['block_size']
        self.file_size = args['file_size']
        self.suffix = args['suffix']
        self.result = args['result']
        self.shellname = args['shellname']
        self.url = args['url']


    def get_block_list(self):
        '''
        根据文件大小和block_size，获取需要切分的block_list
        '''
        size_of_block = float(self.file_size) // float(self.block_size)  # 整除，文件切割区块数 - 1
        block_list = []
        block_list.append((0, self.block_size))

        while len(block_list) < size_of_block:
            start = block_list[-1][-1] + 1
            block_list.append((start, start + self.block_size))

        start = block_list[-1][-1] + 1
        block_list.append((start, ''))

        return block_list


    def get_cmd(self, block_list):
        '''
        生成下载的代码
        '''
        cmd_list = []
        part_list = []  # 用于后期文件合并

        for index,item in enumerate(block_list):
            cmd = 'curl --range {item[0]}-{item[1]} -o {self.suffix}.part{index} {self.url}'.format(**locals())
            part = '{self.suffix}.part{index}'.format(**locals())
            cmd_list.append(cmd)
            part_list.append(part)
        return cmd_list, part_list


    def write_download_shell(self, cmd_list):
        with open(self.shellname, 'w') as fw:
            for cmd in cmd_list:
                fw.write('{}\n'.format(cmd))


    def write_merge_shell(self, part_list):
        with open(self.shellname + '.cat', 'w') as fw:
            fw.write('cat ')
            for part in part_list:
                fw.write('{} '.format(part))
            fw.write(' > {self.result}'.format(**locals()))


    def start(self):
        block_list = self.get_block_list()
        cmd_list, part_list = self.get_cmd(block_list)
        self.write_download_shell(cmd_list)
        self.write_merge_shell(part_list)



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='切割大数据下载')
    parser.add_argument('--block_size', type=int, help='切割文件大小,建议M单位大小')
    parser.add_argument('--file_size', type=float, help='待下载文件大小,误差不要超过一个block_size')
    parser.add_argument('--suffix', help='切割结果前缀')
    parser.add_argument('--result', help='输出文件')
    parser.add_argument('--shellname', help='download的文件名称')
    parser.add_argument('--url', help='带下载文件地址')

    args = vars(parser.parse_args())

    dd = Download(args).start()



