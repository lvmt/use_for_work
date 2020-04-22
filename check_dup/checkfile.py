#!/usr/bin/env python
# -*- coding:utf-8 -*-

import argparse
import re 
from collections import defaultdict


class QC(object):

    def __init__(self, args):
        self.infile = args['infile']
        self.qc_merge_dict = self.get_qc_dict()
        self.items = re.split(";|,|:", args['items'])
        self.values = re.split(";|,|:", args['values'])
        self.check_dict = dict(zip(self.items, self.values))

    def get_qc_dict(self):
        qc_dict = defaultdict(dict)
        with open(self.infile, 'r') as f:
            title_list = f.readline().strip().split('\t')
            for line in f:
                linelist = line.strip().split('\t')
                sam = linelist[0]
                if sam in qc_dict:
                    qc_dict[sam]['rawdata'].append(linelist[4])
                    qc_dict[sam]['effective'].append(linelist[6])
                    qc_dict[sam]['error'].append(linelist[7])
                    qc_dict[sam]['q20'].append(linelist[8])
                    qc_dict[sam]['q30'].append(linelist[9])
                    qc_dict[sam]['gc'].append(linelist[10])
                else:
                    qc_dict[sam]['rawdata'] = [linelist[4]]
                    qc_dict[sam]['effective'] = [linelist[6]]
                    qc_dict[sam]['error'] = [linelist[7]]
                    qc_dict[sam]['q20'] = [linelist[8]]
                    qc_dict[sam]['q30'] = [linelist[9]]
                    qc_dict[sam]['gc'] = [linelist[10]]
        # merge all lane of one sample
        qc_merge_dict = defaultdict(dict)
        for sam in qc_dict:
            rawdata_sum = 0
            effective_sum = 0
            error_sum = 0
            q20_sum = 0
            q30_sum = 0
            gc_sum = 0
            for i in range(len(qc_dict[sam]['gc'])):
                rawdata_sum += float(qc_dict[sam]['rawdata'][i])
                effective_sum += float(qc_dict[sam]['effective'][i]) * float(qc_dict[sam]['rawdata'][i])
                error_sum += float(qc_dict[sam]['error'][i]) * float(qc_dict[sam]['rawdata'][i])
                q20_sum += float(qc_dict[sam]['q20'][i]) * float(qc_dict[sam]['rawdata'][i])
                q30_sum += float(qc_dict[sam]['q30'][i]) * float(qc_dict[sam]['rawdata'][i])
                gc_sum += float(qc_dict[sam]['gc'][i]) * float(qc_dict[sam]['rawdata'][i])
            
            merge_rawdata = rawdata_sum
            merge_effective = round(effective_sum / merge_rawdata, 2)
            merge_error = round(error_sum / merge_rawdata, 2)
            merge_q20 = round(q20_sum / merge_rawdata, 2)
            merge_q30 = round(q30_sum / merge_rawdata, 2)
            merge_gc = round(gc_sum / merge_rawdata, 2)

            qc_merge_dict[sam]['rawdata'] = merge_rawdata
            qc_merge_dict[sam]['effective'] = merge_effective
            qc_merge_dict[sam]['clean'] = round(merge_effective * merge_rawdata / 100, 2)
            qc_merge_dict[sam]['error'] = merge_error
            qc_merge_dict[sam]['q20'] = merge_q20
            qc_merge_dict[sam]['q30'] = merge_q30
            qc_merge_dict[sam]['gc'] = merge_gc

        return qc_merge_dict

    def rawdata(self):
        raw_failed_sams = []
        sams = []
        n = 1
        for sam in self.qc_merge_dict:
            if float(self.qc_merge_dict[sam]['rawdata']) < float(self.check_dict['rawdata']):
                raw_failed_sams.append((n, sam, self.qc_merge_dict[sam]['rawdata']))
                sams.append(sam)
                n += 1
        if raw_failed_sams:
            print("\033[31m \n原始数据小于 {0}G 的样本有 {1}个\n{2}\n\033[0m".format(self.check_dict['rawdata'], n-1, sams))
            print(raw_failed_sams)

    def clean(self):
        clean_failed_sams = []
        n = 1
        sams = []
        for sam in self.qc_merge_dict:
            if float(self.qc_merge_dict[sam]['clean']) < float(self.check_dict['clean']):
                clean_failed_sams.append((n, sam, self.qc_merge_dict[sam]['clean']))
                sams.append(sam)
                n += 1
        if clean_failed_sams:
            print("\033[31m \nclean数据量小于 {0}G 的样本有 {1}个\n{2}\n\033[0m".format(self.check_dict['clean'], n-1, sams))
            print(clean_failed_sams)

    def q20(self):
        q20_failed_sams = []
        n = 1
        sams = []
        for sam in self.qc_merge_dict:
            if float(self.qc_merge_dict[sam]['q20']) < float(self.check_dict['q20']):
                q20_failed_sams.append((n, sam, self.qc_merge_dict[sam]['q20']))
                sams.append(sam)
                n += 1
        if q20_failed_sams:
            print("\033[31m \nq20低于 {0}% 的样本有 {1}个\n{2}\n\033[0m".format(self.check_dict['q20'], n-1, sams))
            print(q20_failed_sams)
    
    def q30(self):
        q30_failed_sams = []
        n = 1
        sams = []
        for sam in self.qc_merge_dict:
            if float(self.qc_merge_dict[sam]['q30']) < float(self.check_dict['q30']):
                q30_failed_sams.append((n, sam, self.qc_merge_dict[sam]['q30']))
                sams.append(sam)
                n += 1
        if q30_failed_sams:
            print("\033[31m \nq30低于 {0}% 的样本有 {1}个\n{2}\n \033[0m".format(self.check_dict['q30'], n-1, sams))
            print(q30_failed_sams)

    def error(self):
        error_failed_sams = []
        n = 1
        sams = []
        for sam in self.qc_merge_dict:
            if float(self.qc_merge_dict[sam]['error']) > float(self.check_dict['error']):
                error_failed_sams.append((n, sam, self.qc_merge_dict[sam]['error']))
                n += 1
                sams.append(sam)
        if error_failed_sams:
            print("\033[31m \nerror高于 {0}% 的样本有 {1}个\n{2}\n \033[0m".format(self.check_dict['error'], n-1, sams))
            print(error_failed_sams)

    def gc(self):
        gc_failed_sams = []
        n = 1
        sams = []
        for sam in self.qc_merge_dict:
            if float(self.qc_merge_dict[sam]['gc']) > float(self.check_dict['gc']):
                gc_failed_sams.append((n, sam, self.qc_merge_dict[sam]['gc']))
                sams.append(sam)
                n += 1
        if gc_failed_sams:
            print("\033[31m \ngc含量高于 {0}% 的样本有 {1}个\n{2}\n \033[0m".format(self.check_dict['gc'], n-1, sams))
            print(gc_failed_sams)

    def check(self):
        check_type = {"rawdata": self.rawdata,
                      "clean": self.clean,
                      "q20": self.q20,
                      "q30": self.q30,
                      "error": self.error,
                      "gc": self.gc
                      }
        for item in self.check_dict:
            check_type[item]()


class Mapping(object):

    def __init__(self, args):
        self.infile = args['infile']
        self.mapping_dict = self.get_mapping_dict()
        self.items = re.split(";|,|:", args['items'])
        self.values = re.split(";|,|:", args['values'])
        self.check_dict = dict(zip(self.items, self.values))

    @staticmethod
    def clean_str(str):
        new_list = [x for x in re.split('\(|\)', str) if x]
        return new_list[-1].strip("%")
    
    def get_mapping_dict(self):
        mapping_dict = defaultdict(dict)
        with open(self.infile, 'r') as f:
            title = f.next()
            try:
                while f:
                    # print(f.next())
                    linelist = f.next().strip().split("\t")
                    mapping_dict[linelist[1]]['dup'] = Mapping.clean_str(linelist[3])
                    mapping_dict[linelist[1]]['mapped'] = Mapping.clean_str(linelist[4])
                    mapping_dict[linelist[1]]['proper_map'] = Mapping.clean_str(linelist[5])
                    mapping_dict[linelist[1]]['pe_map'] = Mapping.clean_str(linelist[6])
                    mapping_dict[linelist[1]]['effective'] = Mapping.clean_str(linelist[12])
                    mapping_dict[linelist[1]]['depth'] = Mapping.clean_str(linelist[14])
                    mapping_dict[linelist[1]]['coverage'] = Mapping.clean_str(linelist[16])
                    mapping_dict[linelist[1]]['20x'] = Mapping.clean_str(linelist[19])
                    mapping_dict[linelist[1]]['10x'] = Mapping.clean_str(linelist[20])
                    mapping_dict[linelist[1]]['4x'] = Mapping.clean_str(linelist[21])
            except StopIteration as e:
                print(e)
        return mapping_dict

    def dup(self):
        dup_failed_sams = []
        n = 1
        for sam in self.mapping_dict:
            if float(self.mapping_dict[sam]['dup']) > float(self.check_dict['dup']):
                dup_failed_sams.append((n, sam, self.mapping_dict[sam]['dup']))
                n += 1
        if dup_failed_sams:
            print("\033[31m\ndup 超过 {}% 的样本有 {} 个 \033[0m".format(self.check_dict['dup'], len(dup_failed_sams)))
            print(dup_failed_sams)

    def mapped(self):
        mapped_failed_sams = []
        n = 1
        for sam in self.mapping_dict:
            if float(self.mapping_dict[sam]['mapped']) < float(self.check_dict['mapped']):
                mapped_failed_sams.append((n, sam, self.mapping_dict[sam]['mapped']))
                n += 1
        if mapped_failed_sams:
            print("\033[31m\nmapped 低于 {}% 的样本有 {} 个\033[0m".format(self.check_dict['mapped'], len(mapped_failed_sams)))
            print(mapped_failed_sams)

    def _4x(self):
        _4x_failed_sams = []
        n = 1
        for sam in self.mapping_dict:
            if float(self.mapping_dict[sam]['4x']) < float(self.check_dict['4x']):
                _4x_failed_sams.append((n, sam, self.mapping_dict[sam]['4x']))
                n += 1
        if _4x_failed_sams:
            print("\033[31m\n4x 低于 {}% 的样本有 {} 个\033[0m".format(self.check_dict['4x'], len(_4x_failed_sams)))
            print(_4x_failed_sams)
        
    def _10x(self):
        _10x_failed_sams = []
        n = 1
        for sam in self.mapping_dict:
            if float(self.mapping_dict[sam]['10x']) < float(self.check_dict['10x']):
                _10x_failed_sams.append((n, sam, self.mapping_dict[sam]['10x']))
                n += 1
        if _10x_failed_sams:
            print("\033[31m\n10x 低于 {}% 的样本有 {} 个\033[0m".format(self.check_dict['10x'], len(_10x_failed_sams)))
            print(_10x_failed_sams)

    def _20x(self):
        _20x_failed_sams = []
        n = 1
        for sam in self.mapping_dict:
            if float(self.mapping_dict[sam]['20x']) < float(self.check_dict['20x']):
                _20x_failed_sams.append((n, sam, self.mapping_dict[sam]['20x']))
                n += 1
        if _20x_failed_sams:
            print("\033[31m\n20x 低于 {}% 的样本有 {} 个\033[0m".format(self.check_dict['20x'], len(_20x_failed_sams)))
            print(_20x_failed_sams)
     
    def effective(self):   #  捕获效率
        eff_failed_sams = []
        n = 1
        for sam in self.mapping_dict:
            if float(self.mapping_dict[sam]['effective']) < float(self.check_dict['effective']):
                eff_failed_sams.append((n, sam, self.mapping_dict[sam]['effective']))
                n += 1
        if eff_failed_sams:
            print("\033[31m\n捕获效率 低于 {}% 的样本有 {} 个\033[0m".format(self.check_dict['effective'], len(eff_failed_sams)))
            print(eff_failed_sams)

    def depth(self):
        depth_failed_sams = []
        n = 1
        for sam in self.mapping_dict:
            if float(self.mapping_dict[sam]['depth']) < float(self.check_dict['depth']):
                depth_failed_sams.append((n, sam, self.mapping_dict[sam]['depth']))
                n += 1
        if depth_failed_sams:
            print("\033[31m\n深度低于 {}% 的样本有 {} 个\033[0m".format(self.check_dict['depth'], len(depth_failed_sams)))
            print(depth_failed_sams)
        
    def coverage(self):
        cov_failed_sams = []
        n = 1
        for sam in self.mapping_dict:
            if float(self.mapping_dict[sam]['coverage']) < float(self.check_dict['coverage']):
                cov_failed_sams.append((n, sam, self.mapping_dict[sam]['coverage']))
                n += 1
        if cov_failed_sams:
            print("\033[31m\n覆盖度低于 {}% 的样本有 {} 个\033[0m".format(self.check_dict['coverage'], len(cov_failed_sams)))
            print(cov_failed_sams)

    def clean_depth(self):
        cdepth_failed_sams = []
        n = 1
        for sam in self.mapping_dict:
            if float(self.mapping_dict[sam]['depth']) - float(self.mapping_dict[sam]['dup']) < float(self.check_dict['cdepth']):
                cdepth_failed_sams.append((n, sam, float(self.mapping_dict[sam]['depth']) - float(self.mapping_dict[sam]['dup'])))
                n += 1
        if cdepth_failed_sams:
            print("\033[31m\n去除dup, 深度低于 {}X的样本有 {}个\033[0m".format(self.check_dict['cdepth'], len(cdepth_failed_sams)))
            print(cdepth_failed_sams)

    def check(self):
        check_type = {"depth": self.depth,
                      "dup": self.dup,
                      "effective": self.effective,
                      "mapped": self.mapped,
                      "4x": self._4x,
                      "10x": self._10x,
                      "20x": self._20x,
                      "coverage": self.coverage,
                      "cdepth": self.clean_depth}
        for item in self.check_dict:
            check_type[item]()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('infile',
                        help="qcstat.xls or mappingstat.xls")
    parser.add_argument('--items',
                        help="which item to check, \
                        ['rawdata', 'clean', 'q20', 'q30', 'error', gc] \
                        ['depth', 'dup', 'effective', 'mapped', '4x', '10x', '20x', 'coverage', 'cdepth']")
    parser.add_argument('--values',
                        help="the limit value")

    args = vars(parser.parse_args())

    if not args['infile']:
        print(parser.print_help())
    elif args['infile'].startswith("qc"):
        qc = QC(args)
        qc.check()
    elif args['infile'].startswith('mapping'):
        map = Mapping(args)
        map.check()
    else:
        print(parser.print_help())

    # qc = QC(args)
    # qc.check()
    # map = Mapping(args)
    # map.check()  

if __name__ == "__main__":
    main()
