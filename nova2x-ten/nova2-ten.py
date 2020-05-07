#!/usr/bin/env python
#-*- coding:utf-8 -*-

""" 
批量进行格式转换
1、sample_list
2、p7接头文件，根据index有所差异
3、p5接头文件，全都是一样的
4、格式转化后的输出路径
"""

import re 
import textwrap
import logging
import colorama
import os

print(__doc__)

class Novo2Xten(object):


    def __init__(self, args):

        self.args = args
        self.sample_list = args['sample_list']
        self.out = args["outdir"]
       

    def logger(self):
        logging.basicConfig(
            format='[%(asctime)s %(funcName)s %(levelname)s %(message)s]',
            datefmt="%Y-%m-%d %H:%M:%S",
            level=logging.INFO)
        log = logging.getLogger(__name__)
        return log


    def print_color(self, text, fg="black", style=""):

        fore_map = {
            'red': colorama.Fore.RED,
            'green': colorama.Fore.GREEN,
            'yellow': colorama.Fore.YELLOW,
            'white': colorama.Fore.WHITE
        }

        style_map = {
            'bright': colorama.Style.BRIGHT,
            'dim': colorama.Style.DIM,
            'normal': colorama.Style.NORMAL
        }

        text = str(text)

        fore = fore_map.get(fg) or fore_map.get("black")
        style = style_map.get(style)
        fore_reset = colorama.Fore.RESET
        style_reset = colorama.Style.RESET_ALL

        text_new = '{fore}{text}{fore_reset}'
        if style:
            text_new = '{style}' + text_new + '{style_reset}'

        return text_new.format(**locals())
      

    def write_p7(self, stat1, pname):
        """ 
        read1对应p7: 需要索引
        read2对应p5：都是一样的,但是还是按照单个样本的index进行写出
        """
        with open(stat1, 'r') as fr, open(pname + '.fa', 'w') as fw:
            for line in fr:
                if line.startswith("ada"):
                    continue
                linelist = line.strip().split("\t")
                p7_str = linelist[3]
                fw.write(p7_str)
                break


    def read_sample_list(self):
        with open(self.sample_list, 'r') as fr:
            for line in fr:
                if line.startswith("#"):
                    continue
                linelist = line.strip().split("\t")
                path = linelist[6]
                flowcell = path.strip().strip('/').split("/")[-1]
                outpath = self.out + "/" + flowcell
                laneid = linelist[0]
                sam = linelist[1]
                libid = linelist[3]
                p7 = re.split(r';|-', linelist[5])[0] 
                p5 = re.split(r';|-', linelist[5])[1]
                
                fq1 = path + "/" + libid + "/" + libid + '_L' + laneid + '_1.fq.gz'             
                fq2 = path + "/" + libid + "/" + libid + '_L' + laneid + '_2.fq.gz'
                adapter1 = path + "/" + libid + "/" + libid + '_L' + laneid + '_1.adapter.list.gz'
                adapter2 = path + "/" + libid + "/" + libid + '_L' + laneid + '_2.adapter.list.gz'
                stat1 = path + "/" + libid + "/" + libid + '_L' + laneid + '_1.adap.stat'
                stat2 = path + "/" + libid + "/" + libid + '_L' + laneid + '_2.adap.stat'

                outfq1 = outpath + "/" + libid + "/" + libid + '_L' + laneid + '_1.fq.gz'
                outfq2 = outpath + "/" + libid + "/" + libid + '_L' + laneid + '_2.fq.gz'
                outad1 = outpath + "/" + libid + "/" + libid + '_L' + laneid + '_1.adapter.list.gz'
                outad2 = outpath + "/" + libid + "/" + libid + '_L' + laneid + '_2.adapter.list.gz'
                
                self.logger().info(self.print_color("开始写出p7文件", 'red', 'dim'))
                self.write_p7(stat1, p7)
                self.logger().info(self.print_color("开始写出p5文件", 'green', 'dim'))
                self.write_p7(stat2, p5)

                cmd = textwrap.dedent("""
                ## read1 进行格式转换
                /TJPROJ4/XJ/script/demultiplex/Bin/fqcheck_adapter_nova2xten_update \\
                    -a {self.out}/{p7}.fa \\
                    -r {fq1} \\
                    -l {adapter1} \\
                    -s {stat1} \\
                    -x {outfq1} \\
                    -d {outad1} \
                
                # read2 进行格式转化
                /TJPROJ4/XJ/script/demultiplex/Bin/fqcheck_adapter_nova2xten_update \\
                    -a {self.out}/{p5}.fa \\
                    -r {fq2} \\
                    -l {adapter2} \\
                    -s {stat2} \\
                    -x {outfq2} \\
                    -d {outad2} \
                """.format(**locals()))
                
                with open(sam + '.sh', 'w') as fw:
                    fw.write(cmd)
                
                self.logger().info(self.print_color(sam + '.sh', 'green'))
                os.system("chmod +x  {shell_name}".format(shell_name = sam + '.sh'))
                self.logger().info(self.print_color("开始生成格式转换脚本", 'yellow'))
                print('')
                


def main():

    demo = Novo2Xten(args)
    demo.read_sample_list()


if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample_list", help="需要进行格式转换的sample_list")
    parser.add_argument("--outdir", help="转换后格式后文件放置位置")
    args = vars(parser.parse_args())
    
    main()

