#!/usr/bin/env python 
#-*- coding:utf-8 -*-

"""
增加case 和 control，
"""


import re 
import logging
from collections import defaultdict
import textwrap
import os



class UPD(object):


    def __init__(self, args):
        self.args = args
        self.sample_info = args['sample_info']
        self.infile = args['infile']  
        self.case = args['case']
        self.control = args['control']
        self.type = args['type']
        self.family_info, self.stat_info = self.get_family_info()  # 先获得每个家系的详细信息


    def logger(self):
        logging.basicConfig(
            format='[%(asctime)s %(funcName)s %(levelname)s %(message)s]',
            datefmt="%Y-%m-%d %H:%M:%S",
            level=logging.INFO)
        log = logging.getLogger(__name__)
        return log


    def safe_open(self, filename, mode):
        try:
            if not filename.endswith(".gz"):
                return open(filename, mode)
            else:
                import gzip
                return gzip.open(filename, mode)
        except IOError:
            print("{filename} does'n exists!!".format(**locals()))


    def get_het(self, genetype):
        genetype_list = re.split(':|-|/', genetype)
        if genetype_list[0] != ".":
            if genetype_list[0] == genetype_list[1] and genetype_list[0] != 0:
                return "hom"
            elif genetype_list[0] == genetype_list[1] and genetype_list[0] == 0:
                return "ref"
            elif genetype_list[0] != genetype_list[1] and "0" in genetype_list:
                return "het"
            else:
                return "chet"
        else:
            return "nocall"  # ./.


    def get_upd_tag(self, gt_proband,gt_father,gt_mother): 
        #gt_proband,gt_father.gt_mother
        """
        根据孩子，父亲，母亲三人的基因型，判断upd类型
        tag : upd类型 + 来自谁
        """
        tag = []  #UN/BI/UA/UI pa/ma
        if self.get_het(gt_proband) in ['ref','hom']:
            if self.get_het(gt_father) in ['ref','hom'] and self.get_het(gt_mother) in ['ref','hom'] and gt_father != gt_mother:
                tag.append('UA')
                if gt_proband == gt_father:
                    tag.append('pa')
                elif gt_proband == gt_mother:
                    tag.append('ma')
            elif (self.get_het(gt_father) in ['ref','hom']) and (gt_proband != gt_father) and (self.get_het(gt_mother) == 'het'):
                tag.append('UI')
                tag.append('ma')
            elif (self.get_het(gt_mother) in ['ref','hom']) and (gt_proband != gt_mother) and (self.get_het(gt_father) == 'het'):
                tag.append('UI')
                tag.append('pa')
            else:
                tag.append('UN')
        elif self.get_het(gt_proband) == 'het':
            if self.get_het(gt_father) in ['ref','hom'] and self.get_het(gt_mother) in ['ref','hom'] and gt_father != gt_mother:
                tag.append('BI')
            else:
                tag.append('UN')
        elif (self.get_het(gt_proband) in ['ref','nocall']) and (self.get_het(gt_father) in ['ref','nocall']) and (self.get_het(gt_mother) in ['ref','nocall']):
            tag.append('noVariant')
        else:
            tag.append('UN')
        #print '_'.join(tag)
        return '_'.join(tag)


    def get_family_info(self):
        """
        得到每个家系的信息, 每个样本在输入文件中的索引
        family_info = {
            'fid':{
                    'proband': {name: name, index: index} 
                    'pa': {name: name, index: index} 
                    'ma': {name: name, index: index}
            },
            'fid':{

            }
        }
        ## 
        stat_info = {
            fid: {
                "1" : {'total':0, 'UN':0, 'UA':0, 'UI':0, 'BI':0},
                "2" : {'total':0, 'UN':0, 'UA':0, 'UI':0, 'BI':0},
            },
            fid:{

            }
        }
        """
        family_info = defaultdict(dict)
        stat_info = defaultdict(dict)
  
        with open(self.sample_info, 'r') as fr:
            for line in fr:
                line = line.strip()
                linelist = line.strip().split("\t")
                nlinelist = [item.lower() for item in linelist]
                if line.startswith('#') and "sampleid" not in nlinelist:
                    continue
                elif "sampleid" in nlinelist:
                    headlist = nlinelist
                    head_index = {head:index for index,head in enumerate(headlist)}
                    self.logger().info('\033[31msample_info表头索引：\n{head_index}\033[0m'.format(**locals())) 
                else:
                    if len(linelist) == len(headlist):
                        famid = linelist[head_index["#familyid"]]
                        probandid = linelist[head_index["sampleid"]]
                        paid = linelist[head_index["pa"]]
                        maid = linelist[head_index["ma"]]

                        family_info[famid]['proband'] = {'name': probandid, 'index': 'None'}
                        family_info[famid]['pa'] = {'name': paid, 'index': 'None'}
                        family_info[famid]['ma'] = {'name': maid, 'index': 'None'}
                    else:
                        self.logger().info("the proband {} has not the information of pa and ma.".format(linelist[head_index['sampleid']]))

        for fid in family_info.keys():
            for chr in range(1, 23):
                stat_info[fid][str(chr)] = {"total": 0, "UN":0, "UA": 0, "UI": 0 , "BI": 0}

        return family_info, stat_info


    def handle_vcf(self):

        """
        从合并的vcf文件中提取每个家系的vcf文件
        """
        case_sams = ",".join([self.family_info[self.case][item]['name'] for item in self.family_info[self.case].keys()])
        control_sams = ",".join([self.family_info[self.control][item]['name'] for item in self.family_info[self.control].keys()])

        cmd = textwrap.dedent('''
        bcftools-1.9 view -s {case_sams} -x -Oz -o {self.case}.{self.type}.vcf.gz {self.infile} 
        bcftools-1.9 view -s {control_sams} -x -Oz -o {self.control}.{self.type}.vcf.gz {self.infile} 
        '''.format(**locals())) 

        self.logger().info('\033[32m处理vcf文件中....\033[0m')
        os.system(cmd)
        self.logger().info('\033[32m处理vcf文件处理完成....\033[0m')


    def analysis_upd(self, fid):
        """
        结果包含2个文件：
        一个原始的文件
        一个统计文件
        """
        fid_infile = '{fid}.{self.type}.vcf.gz'.format(**locals())  # 根据handle_vcf的输出

        with self.safe_open(fid_infile, 'r') as fr, open('.'.join((fid,self.type,'xls')), 'w') as fw, \
            open('.'.join((fid,self.type,'stat.xls')), 'w') as fwstat:
            for line in fr:
                linelist = line.strip().split("\t")
                nlinelist = [x.lower() for x in linelist]
                if line.startswith("##"):
                    continue
                elif 'pos' in nlinelist:
                    headlist = nlinelist
                    format_pos = headlist.index("format")
                    chr_pos = headlist.index("#chrom") if "#chrom" in headlist else headlist.index("chrom")
                    pos_pos = headlist.index("pos")

                    outheader = ["tag", "UPDfam"] + linelist
                    fw.write("\t".join(outheader) + "\n")

                    self.family_info[fid]['proband']['index'] = linelist.index(self.family_info[fid]['proband']['name'])
                    self.family_info[fid]['pa']['index'] = linelist.index(self.family_info[fid]['pa']['name'])
                    self.family_info[fid]['ma']['index'] = linelist.index(self.family_info[fid]['ma']['name'])
                    print('\033[33m家系详细信息\n{}\033[0m'.format(self.family_info[fid]))

                else:
                    gtindex = linelist[format_pos].split(":").index("GT") 
                    tag = []
                    updfam = []
                    chrom = linelist[chr_pos]
                    pos = linelist[pos_pos]

                    if chrom not in ("X", "Y"):
                        gt_proband = linelist[self.family_info[fid]['proband']['index']].split(":")[gtindex]
                        gt_pa = linelist[self.family_info[fid]['pa']['index']].split(":")[gtindex]
                        gt_ma = linelist[self.family_info[fid]['ma']['index']].split(":")[gtindex]
                        upd = self.get_upd_tag(gt_proband, gt_pa, gt_ma)
                        
                        if upd == "noVariant":
                            # print("\033[32m这是个noVariant点，跳过\033[0m")  # 输出 noVariant 点，但是不进行统计
                            continue
                
                        # self.logger().info("该位点的UPD类型为[proband/pa/ma] {upd} {gt_proband} {gt_pa} {gt_ma} || {chrom} {pos} ".format(**locals()))
                        tag.append(upd + '_' + fid)    # upd: UI_famid
                        self.stat_info[fid][chrom]["total"] += 1
                        self.stat_info[fid][chrom][upd.split("_")[0]] += 1    # 类型 +1

                        if upd not in ("UN", "noVariation", "BI"):
                            updfam.append(fid)
                        
                        if updfam != []:
                            outlist = [';'.join(tag), ';'.join(updfam)] + linelist
                        else:
                            outlist = [';'.join(tag), '--'] + linelist
                        fw.write("\t".join(outlist)+"\n")
            
            ## 输出统计文件
            # self.logger().info("====== stat_info 文件加载完成 ======\n{stat_info}".format(**locals()))
            fwstat.write("\t".join(["FamID","CHR","Total","UI","UIperc","UA","UAperc","BI","BIperc","UN","UNperc"])+"\n")
            for ichr in [str(i) for i in range(1, 22+1)]:
                if self.stat_info[fid][ichr]['total'] != 0:
                    uiPerc = self.stat_info[fid][ichr]['UI'] / float(self.stat_info[fid][ichr]['total'])
                    uaPerc = self.stat_info[fid][ichr]['UA'] / float(self.stat_info[fid][ichr]['total'])
                    biPerc = self.stat_info[fid][ichr]['BI'] / float(self.stat_info[fid][ichr]['total'])
                    unPerc = self.stat_info[fid][ichr]['UN'] / float(self.stat_info[fid][ichr]['total'])
                    statout = '%s\t%s\t%d\t%d\t%.4f\t%d\t%.4f\t%d\t%.4f\t%d\t%.4f' %(fid,ichr, self.stat_info[fid][ichr]['total'], \
                        self.stat_info[fid][ichr]['UI'], uiPerc, self.stat_info[fid][ichr]['UA'],  uaPerc, self.stat_info[fid][ichr]['BI'], \
                            biPerc, self.stat_info[fid][ichr]['UN'],unPerc)
                else:
                    statout = '%s\t%s\t0\t0\t-\t0\t-\t0\t-' %(fid,ichr)
                fwstat.write(statout+"\n")
        
        self.logger().info("====== 分析完成 ======")


    def start_analysis(self):

        self.handle_vcf()
        for fid in self.family_info.keys():
            if fid == self.case or fid == self.control:
                self.logger().info('\033[32m开始家系 {fid} 的分析\033[0m'.format(**locals()))
                self.analysis_upd(fid)



def main():
    demo = UPD(args)
    demo.start_analysis()


if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(description="analysis UPD")
    parser.add_argument("--sample_info", help="sample_info文件")
    parser.add_argument("--infile", help="merged vcf文件")
    parser.add_argument('--case', help="case家系ID")
    parser.add_argument('--control', help="control家系ID")
    parser.add_argument('--type', help="输入文件类型，snp还是indel")

    args = vars(parser.parse_args())

    main()

