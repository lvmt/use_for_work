#!/usr/bin/env python 
#-*- coding:utf-8 -*-

import re 
import logging
from collections import defaultdict



class UPD(object):

    def __init__(self, args):
        self.args = args
        self.sample_info = args["sample_info"]
        self.infile = args["infile"]
        self.out = args["out"]

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
            'famid':{
                    'proband': {name: name, index: index} 
                    'pa': {name: name, index: index} 
                    'ma': {name: name, index: index}
                } 
        }

        stat_info = {
            famid: {
                "1" : {'total':0, 'UN':0, 'UA':0, 'UI':0, 'BI':0},
                "2" : {'total':0, 'UN':0, 'UA':0, 'UI':0, 'BI':0},
            }
        }
        """
        family_info = defaultdict(dict)
        stat_info = defaultdict(dict)
  
        self.logger().info(" ========== 开始初始化family_inf和stat_info文件 ==============")

        with open(self.sample_info, 'r') as fr:
            for line in fr:
                line = line.strip()
                linelist = line.strip().split("\t")
                nlinelist = [item.lower() for item in linelist]
                if line.startswith('#') and "sampleid" not in nlinelist:
                    continue
                elif "sampleid" in nlinelist:
                    headlist = nlinelist
                else:
                    if len(linelist) == len(headlist):
                        famid = linelist[headlist.index("#familyid")]
                        probandid = linelist[headlist.index("sampleid")]
                        paid = linelist[headlist.index("pa")]
                        maid = linelist[headlist.index("ma")]

                        family_info[famid]['proband'] = {'name': probandid, 'index': 'None'}
                        family_info[famid]['pa'] = {'name': paid, 'index': 'None'}
                        family_info[famid]['ma'] = {'name': maid, 'index': 'None'}
                    else:
                        self.logger().info("the proband {probandid} has not the information of pa and ma.".format(**locals()))

        for famid in family_info.keys():
            for chr in range(1, 23):
                stat_info[famid][str(chr)] = {"total": 0, "UN":0, "UA": 0, "UI": 0 , "BI": 0}

        self.logger().info('====== family_info 初始化结果 ====== ：\n{family_info}'.format(**locals()))
        self.logger().info('====== stat_info 初始化结果 ====== ：\n{stat_info}'.format(**locals()))
        return family_info, stat_info

    def analysis_upd(self):
        """
        结果包含2个文件：
        一个原始的文件
        一个统计文件
        """
        family_info, stat_info = self.get_family_info()

        with self.safe_open(self.infile, 'r') as fr, open(self.out + '.xls', 'w') as fo, open(self.out + '.stat.xls', 'w') as fostat:
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
                    fo.write("\t".join(outheader) + "\n")

                    for famid in family_info.keys():
                        family_info[famid]['proband']['index'] = linelist.index(family_info[famid]['proband']['name'])
                        family_info[famid]['pa']['index'] = linelist.index(family_info[famid]['pa']['name'])
                        family_info[famid]['ma']['index'] = linelist.index(family_info[famid]['ma']['name'])
                    self.logger().info("======== family_info 字典加载完成，增加样本index信息 =========")
                    print(family_info)

                else:
                    gtindex = linelist[format_pos].split(":").index("GT") 
                    tag = []
                    updfam = []
                    chrom = linelist[chr_pos]
                    pos = linelist[pos_pos]

                    if chrom not in ("X", "Y"):
                        for famid in family_info:
                            gt_proband = linelist[family_info[famid]['proband']['index']].split(":")[gtindex]
                            gt_pa = linelist[family_info[famid]['pa']['index']].split(":")[gtindex]
                            gt_ma = linelist[family_info[famid]['ma']['index']].split(":")[gtindex]
                            upd = self.get_upd_tag(gt_proband, gt_pa, gt_ma)
                            
                            self.logger().info("该位点的UPD类型为[proband/pa/ma] {upd} {gt_proband} {gt_pa} {gt_ma} || {chrom} {pos} ".format(**locals()))
                            tag.append(upd + '_' + famid)    # upd: UI_famid
                            stat_info[famid][chrom]["total"] += 1
                            stat_info[famid][chrom][upd.split("_")[0]] += 1    # 类型 +1

                            if upd not in ("UN", "noVariation", "BI"):
                                updfam.append(famid)
                        
                        if updfam != []:
                            outlist = [';'.join(tag), ';'.join(updfam)] + linelist
                        else:
                            outlist = [';'.join(tag), '--'] + linelist
                        fo.write("\t".join(outlist)+"\n")
            
            ## 输出统计文件
            self.logger().info("====== stat_info 文件加载完成 ======\n{stat_info}".format(**locals()))
            fostat.write("\t".join(["FamID","CHR","Total","UI","UIperc","UA","UAperc","BI","BIperc","UN","UNperc"])+"\n")
            for famid in family_info:
                for ichr in [str(i) for i in range(1, 22+1)]:
                    if stat_info[famid][ichr]['total'] != 0:
                        uiPerc = stat_info[famid][ichr]['UI'] / float(stat_info[famid][ichr]['total'])
                        uaPerc = stat_info[famid][ichr]['UA'] / float(stat_info[famid][ichr]['total'])
                        biPerc = stat_info[famid][ichr]['BI'] / float(stat_info[famid][ichr]['total'])
                        unPerc = stat_info[famid][ichr]['UN'] / float(stat_info[famid][ichr]['total'])
                        statout = '%s\t%s\t%d\t%d\t%.4f\t%d\t%.4f\t%d\t%.4f\t%d\t%.4f' %(famid,ichr, stat_info[famid][ichr]['total'], \
                            stat_info[famid][ichr]['UI'], uiPerc, stat_info[famid][ichr]['UA'],  uaPerc,stat_info[famid][ichr]['BI'], \
                                biPerc,stat_info[famid][ichr]['UN'],unPerc)
                    else:
                        statout = '%s\t%s\t0\t0\t-\t0\t-\t0\t-' %(famid,ichr)
                    fostat.write(statout+"\n")
        
        self.logger().info("====== 分析完成 ======")


def main():
    demo = UPD(args)
    demo.analysis_upd()


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="analysis UPD")
    parser.add_argument("--sample_info", help="sample_info文件，暂时先一个家系一个sample_info吧")
    parser.add_argument("--infile", help="merged vcf文件, 其实其他注释文件也行啊")
    parser.add_argument("--out", help="输出文件前缀，会在原始文件上添加tag，同时生成每条染色体的统计文件")


    args = vars(parser.parse_args())

    main()

