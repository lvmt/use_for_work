#!/usr/bin/env python
#-*- coding:utf-8 -*-

import os
import textwrap
import re
import time


class IGV(object):

    def __init__(self, args):
        self.gene = args.get("gene", "None")
        self.bed = args.get("bed", "None")
        self.projdir = args["projdir"]
        self.sams = args["sams"]
        self.flank = args.get("flank")  # 定义上下游扩展距离，默认2000
        self.outdir = args["outdir"]


    def check_dir(self):
        """
        check dir exists, if not dir , mkdir dir
        """
        if not os.path.exists(self.outdir):
            os.mkdir(self.outdir)


    def gene2bed(self):
        """因为截取文件时，需要输入bed区间，
        如果输入的是gene时，需要得到对应的bed区间
        """
        bed_names = []
        if self.gene:
            bed_names= [x.strip() for x in open(self.gene).readlines()]
            cmd = textwrap.dedent("""
            for i in `less {self.gene}`; do
	        awk -F "\\t" '{{if($1 == "'''$i'''") print $0 }}' /ifs/TJPROJ3/DISEASE/Database/ANNOVAR/humandb/hg19_Genename_from_ncbi.txt | awk -F '\\t' -v OFS='\\t' '{{print$1,$2,$3,$4,$2\":\"$3-{self.flank}\"-\"$4+{self.flank},$2\"_\"$3\"-\"$4}}'>$i.bed
            done""".format(**locals())) 
        os.system(cmd)

        print('\033[32m输入了gene文件, 参与分析的基因名称: {bed_names}\033[0m'.format(**locals()))
        return bed_names


    def bed2bed(self):
        """
        无论是gene还是bed区间，我们最终需要一个输入的bed文件,
        参数实例：4p16:4:966959-3954926
        """
        with open(self.bed, 'r') as fr:
            bed_names = []
            for line in fr:
                line = line.strip()
                name,chr,start,end = re.split(r':|-|_', line)
                newstart = int(start) - self.flank
                newend = int(end) + self.flank
                cmd = textwrap.dedent("""
                echo -e '{name}\t{chr}\t{start}\t{end}\t{chr}:{newstart}-{newend}\t{chr}_{start}-{end}' > {name}.bed
                """.format(**locals())) 

                bed_names.append(name)
                os.system(cmd)
        print('\033[32m输入了bed文件, 参与分析的基因名称: {bed_names}\033[0m'.format(**locals()))
        return bed_names


    def check_bed(self, bedfile):
        """有时候根据gene提取bed文件时，无法提取，
        生成空的bed文件，但是程序还是会接着运行，很明显这是有问题的
        """
        if os.path.exists(bedfile):
            if os.path.getsize(bedfile):
                return True 
            else:
                print("{bedfile}是空的, 请检查".format(**locals()))
                return False
        else:
            print('{bedfile}文件不存在, 请check...'.format(**locals()))
            return False


    def start(self):
        self.check_dir()
        bednames = []
        if self.gene:
            bednames += self.gene2bed()
        if self.bed:
            bednames += self.bed2bed() 
        bednames = " ".join(bednames)

        print('\033[33m最终参与分析的基因名称: {bednames}\033[0m'.format(**locals()))
        print('\033[33m最终参与分析的样本名称：{0}\033[0m'.format([x.strip() for x in open(self.sams).readlines()]))
        time.sleep(2)

        # check bed文件有效性
        not_exist = []
        exist = []
        for name in bednames.split(" "):
            bedfile = name + ".bed"
            if self.check_bed(bedfile):
                exist.append(name)
            else:
                not_exist.append(name)

        print("\033[32mbed文件有效性check中.....\033[0m")
        if not not_exist:
            print("\033[33mbed文件全部有效: {exist}\033[0m".format(**locals()))
        else:
            print("\033[33m有效bed文件: {exist}\033[0m".format(**locals()))
            exit("\033[363无效bed文件：{not_exist}\033[0m".format(**locals()))

        cmd = textwrap.dedent("""
        for i in `less {self.sams}`;do
	        for j in {bednames};do
                pos=`cut -f5 $j.bed`
                # name=`cut -f6 $j.bed`
                samtools view -b -S -h {self.projdir}/Mapping/$i.$i/$i.final.bam $pos >$i.$j.bam;
                samtools index $i.$j.bam;
                /ifs/TJPROJ3/DISEASE/share/Software/Circos/IGVTools/igvtools count -w 10 --query $pos {self.projdir}/Mapping/$i.$i/$i.final.bam $i.$j.wig novo37;
                #/ifs/TJPROJ3/DISEASE/share/Software/Circos/IGVTools/igvtools toTDF $i.$name.wig $i.$name.tdf novo37;
	        done
        done
        """.format(**locals())) 

        os.system(cmd)


if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(description="\033[32m对付广妇儿截图项目\033[0m")
    
    parser.add_argument('-p', '--projdir', help="分析项目路径")
    parser.add_argument('-s', '--sams', help="需要截取bam的名称,文件")
    parser.add_argument('-o', '--outdir', help="结果文件输出路径", default=os.getcwd())
    parser.add_argument('-g', '--gene', help="需要截取的gene名称,文件")
    parser.add_argument('-b', '--bed', help="需要截取的bed区间,文件|实例:4p23:4:start:end")
    parser.add_argument('-f', '--flank', help="截取上下游扩展距离, default=%(default)s", default=2000)
    
    args = vars(parser.parse_args())

    print(args)

    igv = IGV(args)
    igv.start()
