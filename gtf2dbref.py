#!/usr/bin/env python3
#-*- coding:utf-8 -*-


'''
将gff文件转换为dbref格式
'''
import re
import os


class SimpleFile:
    # 简化gff文件
    def __init__(self, gff, outfile):
        self.gff = gff
        self.outfile = outfile


    def start(self):
        with open(self.gff, 'r') as fr, open(self.outfile, 'w') as fw:
            flag_tran = ''
            flag = region = ''  # 标签

            for line in fr :
                linelist = line.strip('\n').split('\t')
                region = linelist[2]  # gene,mRNA,exon,CDS
                if region in ('gene', 'mRNA'):
                    continue
                _chr = linelist[0]
                if _chr.startswith('NT'):
                    continue
                _chr = re.search(r'NC_0*([1-9]\d*)\.', _chr).group(1)
                _chr = 'chr' + _chr
                start = linelist[3]
                end = linelist[4]
                strand = linelist[6]
                info = linelist[8]
                gene = re.search(r'gene=(.*?);', info).group(1)
                transcript = re.search(r'Parent=rna-(.*?);', info).group(1)
                region = linelist[2]  # exon CDS

                # core code
                if (not transcript == flag_tran) or (not flag_region == region):
                    flag_tran = transcript
                    flag_region = region
                    num = 0

                num += 1

                fw.write('{_chr}\t{start}\t{end}\t{gene}\t{strand}\t{transcript}\t{region}-{num}\n'.format(**locals()))



class GetRelationShip:
    # 将exon和CDS关系进行映射

    def __init__(self, infile, outfile):
        self.infile = infile
        self.outfile = outfile


    def split_region(self, exon_start, exon_end, cds_start, cds_end, strand, exon, cds):
            if exon_start == cds_start and exon_end == cds_end:
                return [(exon_start, exon_end, exon, cds)]
            elif exon_end == cds_end and strand == '+':  # 5'UTR
                region1 = (exon_start, int(cds_start) - 1, exon, "5'UTR")
                region2 = (cds_start, cds_end, exon, cds)
                return [region1, region2]
            elif exon_start == cds_start and strand == '+':  # 3'UTR
                region1 = (cds_start, cds_end, exon, cds)
                region2 = (int(cds_end) + 1, exon_end, exon, "3'UTR")
                return [region1, region2]
            elif exon_start == cds_start and strand == '-':  # 5'UTR
                region1 = (cds_start, cds_end, exon, cds)
                region2 = (int(cds_end) + 1, exon_end, exon, "5'UTR")
                return [region2, region1]
            elif exon_end == cds_end and strand == '-':  # 3'utr
                region1 = (exon_start, int(cds_start) - 1, exon, "3'UTR")
                region2 = (cds_start, cds_end, exon, cds)
                return [region2, region1]
            else:
                # 最后一个外显子，上面既有5UTR-exon-3UTR
                if strand == '-':
                    region1 = (exon_start, int(cds_start) - 1, exon, "3'UTR")
                    region2 = (cds_start, cds_end, exon, cds)
                    region3 = (int(cds_end) + 1, exon_end, exon, "5'UTR")
                    return [region3, region2, region1]
                else:
                    region1 = (exon_start, int(cds_start) - 1, exon, "5'UTR")
                    region2 = (cds_start, cds_end, exon, cds)
                    region3 = (int(cds_end) + 1, exon_end, exon, "3'UTR")
                    return [region1, region2, region3]


    def start(self):
        with open(self.infile, 'r') as fr, open(self.outfile, 'w') as fw:
            for line in fr:
                linelist = line.strip('\n').split('\t')
                _chr = linelist[0]
                exon_start = linelist[1]
                exon_end = linelist[2]
                gene = linelist[3]
                strand = linelist[4]
                transcript = linelist[5]
                exon = linelist[6]  # 为.的话，表明UTR区域
                cds_start = linelist[8]
                cds_end = linelist[9]
                cds = linelist[13]

                if cds == '.':
                    fw.write('{_chr}\t{exon_start}\t{exon_end}\t{gene}\t{transcript}\t{strand}\t{exon}\tUTR\n'.format(**locals()))
                else:
                    new_region = self.split_region(exon_start, exon_end, cds_start, cds_end, strand, exon, cds)
                    for (start, end, exon, cds) in new_region:
                        fw.write('{_chr}\t{start}\t{end}\t{gene}\t{transcript}\t{strand}\t{exon}\t{cds}\n'.format(**locals()))



class Gtf2bdref:
    # 转换主程序
    # 中间会生成一些中间文件,建议在独立目录下运行程序
    # 最终决定单个转录本进行计算，是由于会出现意外情况(NM_001282672.2)
    # 由于是单个转录本的运行，导致程序运行时间较长，但是结果的准确性达到100%。
    def __init__(self, gff, dbref):
        self.gff = gff
        self.dbref = dbref


    def start(self):
        simple_gff = self.gff + '.simple'
        SimpleFile(self.gff, simple_gff).start()
        # 获取全部转录本文件
        os.system('less {simple_gff} | cut -f 6 | sort | uniq > all.transcripts'.format(**locals()))
        with open('all.transcripts', 'r') as fr:
            for line in fr:
                tran = line.strip('\n')
                os.system('grep {tran} {simple_gff} | grep exon > tmp.exon'.format(**locals()))
                os.system('grep {tran} {simple_gff} | grep CDS > tmp.cds'.format(**locals()))
                os.system('bedtools intersect -a tmp.exon -b tmp.cds -wao > tmp.exon.cds')

                GetRelationShip('tmp.exon.cds', 'tmp.dbref').start()
                os.system('cat tmp.dbref >> {self.dbref}'.format(**locals()))



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='\033[1;32mTransvar gff to dbref\033[0m')
    parser.add_argument('--gff', help='input file, gff format')
    parser.add_argument('--dbref', help='outfile, dbref format')

    args = parser.parse_args()

    Gtf2bdref(args.gff, args.dbref).start()


