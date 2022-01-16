#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   command.py
@Time    :   2021/11/13 11:40:43
@Author  :   lvmt 
@Version :   1.0
'''

# here put the import lib
'''
模块主程序
'''


from simplify.simplify_vep_annotation import SimpleVep
from AddDatabase.add_database import DataBaseAnnotation
from Filter.filter import FilterPipeline
from logicalcheck.logical_check import LogicalCheck
        
        
        
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description='\033[1;32mBGI注释过滤模块\ncontact: lvmengting@genomics.cn\033[0m'
    )
    subparser = parser.add_subparsers(
        title = 'sub-command',
        dest='subparser_name',
        metavar=''
    )
    
    
    #vep简化 子命令
    simplify = SimpleVep()
    subparser_simplify = subparser.add_parser(
        'simplify',
        formatter_class=argparse.RawTextHelpFormatter,
        help='VEP annotation simplify'
    )
    simplify.parser_vep_simplify(subparser_simplify)


    #额外数据库注释 子命令
    database = DataBaseAnnotation()
    subparser_database = subparser.add_parser(
        'database_annotation', 
        formatter_class=argparse.RawTextHelpFormatter,
        help='额外数据库注释,(cosmic, uniport..)'
    )
    database.parser_database_annotation(subparser_database)
    
    
    #过滤pipeline 子命令
    filter = FilterPipeline()
    subparser_filter = subparser.add_parser(
        'filter',
        formatter_class=argparse.RawTextHelpFormatter,
        help='过滤步骤'
    )
    filter.parser_filter(subparser_filter)
    
    
    #其他人工逻辑判断
    check = LogicalCheck()
    subparser_check = subparser.add_parser(
        'check',
        formatter_class=argparse.RawTextHelpFormatter,
        help='逻辑判断'
    )
    check.parser_check(subparser_check)
    
    

    args = parser.parse_args()
    print(args)
    args.func(**vars(args))    
    
        
        
