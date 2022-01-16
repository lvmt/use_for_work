#!/usr/bin/env python
#-*- coding:utf-8 -*-
'''
@Author: lvmengting
@Date: 2021-11-14 13:40:30
@Last Modified by:   lvmengting  
@Last Modified time: 2021-11-14 13:40:30
'''


from Filter.filter_func import FilterFunc
from Filter.filter_reads import FilterReads
from Filter.filter_localcontrol import FilterLocalControl
from Filter.filter_database import FilterDataBase
from Filter.filter_synonymy import FilterSynonymy



class FilterPipeline:
    '''
    对注释结果做系列过滤
    '''
    def parser_filter(self, parser):
        parser.add_argument(
            '--infile', help='输入文件'
        )
        parser.add_argument(
            '--result_suffix', help='结果文件前缀'
        )

        func_group = parser.add_argument_group(title='功能过滤')
        func_group.add_argument(
            '--func_config', help='功能过滤的配置文件'
        )

        database_group = parser.add_argument_group(title='本地数据库过滤')
        database_group.add_argument(
            '--database_config', help='本地数据库配置文件'
        )

        parser.add_argument(
            '--filter_synonymy', action='store_true', help='进行同义突变过滤'
        )
        
        read_group = parser.add_argument_group(title='基于read的过滤')
        read_group.add_argument(
            '--reads', type=int, help='read支持数过滤'
        )

        local_group = parser.add_argument_group(title='基于本地freq数据库进行过滤')
        local_group.add_argument(
            '--local_snv', help='本地snv频率数据库'
        )
        local_group.add_argument(
            '--local_indel', help='本地indel频率数据库'
        )

        parser.set_defaults(func=self.start)

        # freq_group = parser.add_argument_group(title='人群频率数据库过滤')
        # freq_group.add_argument(
        #     '--freq', type=float, help='频率阈值',
        # )
        # freq_group.add_argument(
        #     '--freq_database_ratio', type=float,
        #     help='至少多少个频率库小于阈值'
        # )

        # score_group = parser.add_argument_group(title='基于分值的过滤')
        # score_group.add_argument(
        #     '--score', help='分数过滤'
        # )
        
        
    def start(self, **args):

        if args['func_config']:
            FilterFunc(**args).start()
        
        if args['database_config']:
            FilterDataBase(args).start()
                
        if args['reads']:
            FilterReads(args).start()

        # if args['freq']:
        #     FilterFreqDatabase(**args).start()      

        if args['filter_synonymy']:
            FilterSynonymy(args).start()

        if all(
            (args.get('infile'),
            args.get('result_suffix'),
            args.get('local_snv'),
            args.get('local_indel'))
            ):
            FilterLocalControl(args).start()


        print('this is start function')
        print(args)

        