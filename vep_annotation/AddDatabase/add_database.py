#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   add_database.py
@Time    :   2021/11/13 11:40:43
@Author  :   lvmt 
@Version :   1.0
'''


from collections import defaultdict
from AddDatabase.uniport import UNIPORT
from AddDatabase.cosmic import COSMIC
from AddDatabase.target_gene import TargetGene
from AddDatabase.maploc import MapLoc
from AddDatabase.tmb import TMB
from AddDatabase.func_region import FuncRegion



class DataBaseAnnotation:
    '''
    给注释结果增加其他的数据库注释字段
    '''
    def parser_database_annotation(self, parser):
        parser.add_argument(
            '--infile', help='输入文件'
        )
        parser.add_argument(
            '--result', help='输出文件'
        )
        uniport_group = parser.add_argument_group(title='增加uniport数据库注释')
        uniport_group.add_argument(
            '--uniport', help='uniport库文件'
        )

        cosmic_group = parser.add_argument_group(title='增加cosmic数据库注释')
        cosmic_group.add_argument(
            '--cosmic', help='cosmic库文件'
        )

        target_gene_group = parser.add_argument_group(title='增加target_gene字段注释')
        target_gene_group.add_argument(
            '--target_gene_pos', help='包含chgvs和phgvs库文件'
        )
        target_gene_group.add_argument(
            '--target_gene_func', help='包含func信息库文件'
        )
        target_gene_group.add_argument(
            '--target_gene_exon', help='包含特定基因的库文件'
        )

        tmb_group = parser.add_argument_group(title='增加TMB字段')
        tmb_group.add_argument(
            '--driver_list', help='驱动基因列表'
        )
        tmb_group.add_argument(
            '--driver_yaml', help='特例基因配置文件',
        )
        tmb_group.add_argument(
            '--tmb_freq', help='tmb突变阈值', type=float, default=1.5
        )
        tmb_group.add_argument(
            '--chip_size', help='暂时用不上'
        )

        parser.add_argument(
            '--maploc', help='变异条带区注释'
        )

        parser.add_argument(
            '--func_region_relation', help='cds和exon对应关系配置文件'
        )
        
        parser.set_defaults(func=self.start)
        
     
    def start(self, **args):
        if args['uniport']:
            UNIPORT(args).start()

        if args['cosmic']:
            COSMIC(args).start()

        if args['target_gene_pos'] and args['target_gene_func'] and args['target_gene_exon']:
            TargetGene(args).start()

        if args['maploc']:
            MapLoc(args).start()

        if all((args.get('infile'),
               args.get('result'),
               args.get('driver_list'),
               args.get('driver_yaml'),
               args.get('tmb_freq')
               )):
               
            TMB(args).start()

        if args['func_region_relation']:
            FuncRegion(args).start()


        







    
    
    
    
    