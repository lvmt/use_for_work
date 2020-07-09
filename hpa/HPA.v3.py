#!/usr/bin/env python
#-*- coding:utf-8 -*- 

import re
import os
from collections import OrderedDict,defaultdict


Usage = """\033[1;36m
单倍型分析：
单倍型分析，只能针对SNP位点
1、指定区域内的SNP位点单倍型分析
2、指定rs号的单倍型分析（某个SNP 和 某个SNP），必须是同一条染色体上位点
3、指定基因的单倍型分析（分析某个基因的SNP位点是否符合单倍型分析）
4、指定chr:pos的单倍型分析

分析模式：
根据样本的特性：可以是散发单倍型分析，也可以是家系样本单倍型分析
-assocCC: 散发样本
-assocTDT: 家系样本

输入文件： 2个
snp.vcf文件
snp.vcf注释文件

图片展示问题
如果图片展示不全，可能是进行单倍型分析中的info文件的第一列名字太长了
\033[0m"""


class HPA(object):

    def __init__(self, args):
        self.args = args
        self.vcf = args['vcf']
        self.anno = args['anno']
        self.gene = args['gene']
        self.region = args['region']
        self.rs = args['rs']
        self.site = args['site']
        self.info = args['info']
        self.out = args['out']
        self.analysis = args['analysistype']

    def get_head_index(self, headlist):
        '''
        得到头文件每列对应的索引值
        '''
        index_dict = OrderedDict
        index_dict = {value.lower():index for index,value in enumerate(headlist)}
        return index_dict
    
    def safe_open(self, file, mode='r'):
        if not file.endswith('.gz'):
            return open(file)
        elif file.endswith('.gz'):
            import gzip
            return gzip.open(file)
        else:
            exit('wrong file type')
            
    def get_need_vcf(self):
        '''
        处理gene模式|处理region模式|处理指定rs模式,
        从注释文件中提取目的位置信息及rs号，
        然后在去vcf文件里面提取，并增加对应的rs号信息
        anno_dict = {
            '1_123': rs,
            '1_345': rs1
        }
        '''
        anno_dict = {}
        with self.safe_open(self.anno, 'r') as fr:
           for line in fr:
                linelist = line.strip().split('\t')
                if linelist[0].lower() in ('chrom', 'priority'):
                    headlist = linelist
                    continue
                index_dict = self.get_head_index(headlist)
                gene_index = index_dict['genename']
                chr_index = index_dict['chrom'] or index_dict['chr']
                pos_index = index_dict['pos']
                rs_index = index_dict['id']

                k = linelist[chr_index] + '_' + linelist[pos_index]
                v = linelist[rs_index]

                if self.gene:
                    # print('\033[32m进行指定基因的HPA分析\033[0m')
                    if self.gene == linelist[gene_index]:
                        anno_dict[k] = v
                    else:
                        pass

                if self.region: # chr:start:end
                    # print('\033[32m进行指定区间的HPA分析\033[0m')
                    chr,start,end = re.split(r':|,|-', self.region)
                    if chr == linelist[chr_index]:
                        if int(linelist[pos_index]) >= int(start) and int(linelist[pos_index]) <= int(end):
                            anno_dict[k] = v
                
                if self.rs:
                    # print('\033[32m进行指定rs号的HPA分析')
                    all_rs = re.split(r':|,|-|_|\||', self.rs)
                    if linelist[rs_index] in all_rs:
                        anno_dict[k] = v

                if self.site:
                    # 指定染色体和pos， chr:pos:pos2:pos3
                    all_site = re.split(r':|_|,|-', self.site)
                    new_all_site = [ "_".join((all_site[0], item)) for item in all_site[1:]]
                    if k in new_all_site:
                        anno_dict[k] = v

        print('\033[32m打印anno_dict内容\033[0m')
        print(anno_dict)
        
        with self.safe_open(self.vcf, 'r') as fr, open(self.out+'.vcf', 'w') as fw:
            for line in fr:
                linelist = line.strip().split('\t')
                if line.startswith('##'):
                    fw.write(line)
                elif linelist[0].lower() in ('#chrom', '#chr'):
                    fw.write(line)
                    headlist = linelist
                    index_dict = self.get_head_index(headlist)
                    chr_index = index_dict['#chrom']
                    pos_index = index_dict['pos']
                    rs_index = index_dict['id']
                    ref_index = index_dict['ref']
                    alt_index = index_dict['alt']
                else:
                    k = linelist[chr_index] + '_' + linelist[pos_index]
                    if k in anno_dict:
                        if anno_dict[k] != '.':
                            linelist[rs_index] = anno_dict[k]
                        else:
                            linelist[rs_index] = "_".join((linelist[chr_index], linelist[pos_index],\
                                                           linelist[ref_index], linelist[alt_index]))  
                        fw.write("\t".join(linelist) + '\n') 
                    else:
                        pass
 
    def get_ped_use_plink(self):
        '''单倍型分析中需要ped文件：
        ped文件包含2部分内容：
        1、利用plink从vcf文件中获取原始的ped文件
        2、从输入的样本信息文件(包含型别，患病信息等)，更新1中ped文件的前6列内容，
        如果不存在信息文件，1作为最终结果
        如果存在信息文件，2作为最终结果
        ped 文件格式：
        famid   individualid    paid    maid    sex phenotype
        '''
        # 利用vcf文件得到ped文件
        vcf2ped_cmd = 'plink --vcf {out}.vcf --recode --make-bed --out {out}.vcf.plink'.format(**args)
        os.system(vcf2ped_cmd)
        info_dict = self.get_info_dict() if self.info else {}
        if info_dict:
            with open(self.out + '.vcf.plink.ped', 'r') as fr, open(self.out + '.vcf.final.plink.ped', 'w') as fw:
                for line in fr:
                    linelist = line.strip().split(" ")
                    samid = linelist[1]
                    if samid in info_dict:
                        linelist[0] = info_dict[samid]['famid']
                        linelist[2] = info_dict[samid]['pa']
                        linelist[3] = info_dict[samid]['ma']
                        linelist[4] = info_dict[samid]['sex']
                        linelist[5] = info_dict[samid]['phenotype']
                    else:
                        print('\033[32样本不在输入的info文件中，请检查\033[0m')
                    fw.write("\t".join(linelist) + '\n')
        else:
            with open(self.out + '.vcf.plink.ped', 'r') as fr, open(self.out + '.vcf.final.plink.ped', 'w') as fw:
                for line in fr:
                    linelist = line.strip().split(' ')
                    linelist[5] = '0'   # 默认输出的pheno为-9，后续无法识别出来
                    fw.write('\t'.join(linelist) + '\n')

    def change_sex(self, sex):
        if sex.lower() in ('m', 'male'):
            return '1'
        if sex.lower() in ('f', 'female'):
            return '2'
        if sex.lower() in ('u', 'unknown'):
            return '3'

    def change_pheno(self, pheno):
        if pheno.lower() in ('n'):
            return '1'
        if pheno.lower() in ('p'):
            return '2'
        if pheno.lower() in ('u', 'unknown'):
            return '0'

    def get_info_dict(self):
        '''从外部输入信息文件中，得到样本的父母信息，型别信息，患病信息等
        然后对get_ped的初始ped文件进行更新
        '''
        info_dict = defaultdict(dict)
        with open(self.info, 'r') as fr:
            for line in fr:
                linelist = line.strip().split('\t')
                if linelist[2].lower() in ('sex', 'SEX'):
                    continue
                famid = linelist[0]
                samid = linelist[1]
                sex = self.change_sex(linelist[2])
                phenotype = self.change_pheno(linelist[3])
                pa = linelist[4] if len(linelist) > 4 else '0'
                ma = linelist[5] if len(linelist) > 4 else '0'

                info_dict[samid] = {
                    'famid': famid,
                    'sex': sex,
                    'phenotype': phenotype,
                    'pa': pa,
                    'ma': ma
                }
        #print(info_dict)
        return info_dict

    def bim2info(self):
        '''从plink结果的bim文件中，提取指定行，得到单倍型分析中需要的info文件
        '''
        with open(self.out + '.vcf.plink.bim', 'r') as fr, open(self.out + '.vcf.HWE.info', 'w') as fw:
            for line in fr:
                linelist = line.strip().split('\t')
                fw.write("\t".join((linelist[1], linelist[3])) + '\n')

    def hap(self):
        '''
        单倍型分析：
        单倍型分析有很多可选参数，记得添加灵活接口参数
        '''
        pedfile = self.out + '.vcf.final.plink.ped'
        info = self.out + '.vcf.HWE.info'
        result = self.out + '.result'
        analysis = self.analysis
        cmd = '''
        java -jar /WORK/Disease/yincan/packages/Haploview.jar -pedfile {pedfile} -info {info} \
        -memory 20000 -maxdistance 14000 -png -out {result} -dprime -nogui -minMAF 0.05 -ldvalues RSQ {analysis}
        '''.format(**locals())
        os.system(cmd)

    def start(self):
        self.get_need_vcf()
        self.get_ped_use_plink()
        self.bim2info()
        self.hap()



if __name__ == '__main__':
    
    import argparse
    
    parser = argparse.ArgumentParser(description=Usage, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-r', '--region', help='进行HPA分析的区间')
    parser.add_argument('-g', '--gene', help='进行HPA分析的gene')
    parser.add_argument('-s', '--site', help='根据chr+pos进行HPA分析:chr:pos1:pos2...')
    parser.add_argument('--rs', help='对指定rs进行HPA分析，需要在同一条染色体上')
    parser.add_argument('-v', '--vcf', help='输入的vcf文件')
    parser.add_argument('-a', '--anno', help='输入为vcf的注释文件')
    parser.add_argument('-f', '--info', help='包含样本型别，患病信息的配置文件')
    parser.add_argument('-o', '--out', help='输出文件前缀', default='result')
    parser.add_argument('-t', '--analysistype', help='散发还是家系样本', choices=['-assocCC', '-assocTDT'], default='-assocCC')
    
    args = vars(parser.parse_args())
    
    print(args)
    
    hpa = HPA(args)
    # hpa.get_need_vcf()
    # hpa.get_ped_use_plink()
    # hpa.get_info_dict()
    # hpa.get_ped_use_plink()
    # hpa.bim2info()
    hpa.start()
    
