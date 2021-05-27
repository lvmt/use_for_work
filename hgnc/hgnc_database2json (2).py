#!/usr/bin/python3
# -*- coding:utf-8 -*-
# @Author lvmengting
# @Time   2021/5/27 11:28
# @Email  13554221497@163.com
# @File   hgnc_database2json.py



"""
将hgnc的数据库文件转换为json,存储本地,减少查询时间
"""

from collections import defaultdict
import pandas as pd
import numpy as np
import json


class Database2Json:

    def __init__(self, args):
        self.database = args['database']
        self.jsonfile = args['jsonfile']

    def get_database_info(self):
        """
        利用dataframe读取数据
        approved symbol：只有一个gene名
        previous symbols：可能存在多个基因名
        :return:
        """
        database_info = {}
        df = pd.read_csv(self.database, sep='\t')
        df = df[df['Status'] != 'Symbol Withdrawn']  # 剔除有信心的标准名的基因
        col = ['Approved symbol', 'Previous symbols', 'Status', 'Date approved']  # 提取指定列
        df = df[col]
        df = df.fillna('none')

        for item in df.values:
            symbol = item[0]
            previous = item[1]
            status = item[2]
            date = item[3]
            # print(type(previous), 'hhhhhh')
            genes = [symbol] + previous.split(', ') if not previous == 'none' else [symbol]  # 注意pandas里面空值的类型.
            for gene in genes:
                database_info[gene] = {
                    'symbol': symbol,
                    'previous': previous,
                    'status': status,
                    'date': date
                }
        return database_info


    def info2json(self, database_info):
        """
        将数据存储值本地json文件
        :return:
        """
        with open(self.jsonfile, 'w') as fw:
            fw.write(json.dumps(database_info, ensure_ascii=False, indent=4))


    def start(self):
        database_info = self.get_database_info()
        self.info2json(database_info)



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='change database txt 2 local json')
    parser.add_argument('--database', help='hgnc database.txt file')
    parser.add_argument('--jsonfile', help='output json file')

    args = vars(parser.parse_args())
    dd = Database2Json(args)
    dd.start()
