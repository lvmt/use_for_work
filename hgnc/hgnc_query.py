#!/usr/bin/python3
# -*- coding:utf-8 -*-
# @Author lvmengting
# @Time   2021/5/27 14:35
# @Email  13554221497@163.com
# @File   hgnc_query.py


"""
查询基因的标准名
"""


import re
import os
import json
import textwrap

class QueryHgnc:

    def __init__(self, args):
        self.query = args['query']
        self.json = args['json']
        self.show = args['show']
        self.result = args['result']

    def get_database_info(self):
        """
        将json文件转换为dict_info文件
        :return:
        """
        with open(self.json, 'r', encoding='utf-8') as fr:
            database_info = json.load(fr)
        return database_info

    def get_query_list(self):
        """
        获取详细的查询基因列表
        :return:
        """
        query_list = []
        if os.path.exists(self.query):
            with open(self.query, 'r') as fr:
                for gene in fr:
                    query_list.append(gene.strip())
        else:
            query_list = re.split(r',|;', self.query)
        return query_list

    def get_query_result(self, databse_info, query_list):
        """
        获取全部的查询结果
        :param databse_info:
        :param query_list:
        :return:
        """
        query_result = []
        invalid_gene = []  # 没有记录在HGNC库文件中,应该是书写错误
        for query in query_list:
            result = databse_info.get(query)
            if result:
                result['query'] = query
                query_result.append(result)
            else:
                invalid_gene.append(query)
        if invalid_gene:
            print('\033;1;32m无效基因\n{}'.format('|'.join(invalid_gene)))
            exit()
        return query_result


        pass

    def show_query_result(self, query_result):
        """
        打印查询结果
        :param query_result:
        :return:
        """
        for result in query_result:
            tmp = textwrap.dedent('''
            \033[1;32mquery_gene: {query}
            \033[1;33mhgnc_gene : {symbol}
            \033[1;33mprevious  : {previous}
            \033[1;33mstatus    : {status}
            \033[1;33mdate      : {date}
            '''.format(**result))
            print(tmp)

    def save_query_result(self, query_result):
        """
        将查询结果输出到excel
        :param query_result:
        :return:
        """
        with open(self.result, 'w') as fw:
            fw.write('query\tHGNC_gene\tprevious_symbol\tstatus\tdate\n')
            for result in query_result:
                tmp = '{query}\t{symbol}\t{previous}\t{status}\t{date}\n'.format(**result)
                fw.write(tmp)

    def start(self):
        database_info = self.get_database_info()
        query_list = self.get_query_list()
        query_result = self.get_query_result(database_info, query_list)
        if self.show:
            self.show_query_result(query_result)
        self.save_query_result(query_result)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Query HGNC about Gene')
    parser.add_argument('--query', help='query gene list of file.')
    parser.add_argument('--json', help='HGNC database json file')
    parser.add_argument('--result', help='output query file', default='hgnc_query_result')
    parser.add_argument(
        '--show',
        help='if print query result',
        action='store_true')

    args = vars(parser.parse_args())
    qq = QueryHgnc(args)
    qq.start()
