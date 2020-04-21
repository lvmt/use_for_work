#!/usr/bin/env python3
#-*- coding:utf-8 -*-

"""检测ONT下机数据
"""

import sys 
import os
import time
sys.path.append('/ifs/TJPROJ3/DISEASE/share/Software/python/site-packages')
import colorama



"""需要的文件
1、类似sample_list 的文件
##lib  flowcell  path name
2、项目路径

"""


## get date ; 区分脚本执行日期
suffix = time.strftime("%Y_%m_%d", time.localtime())


def add_color(text, color=colorama.Fore.CYAN):
    return '%s%s%s' % (color, text, colorama.Fore.RESET)

# def checkfile(query):
#     if os.path.exists(query):
#         return True
#     return False

## count number
done = 1
undone = 1
wrong = 1
    
def  writeshell(samplelist, projdir):
    """samplelist
    
    lib cell path name
    """
    donelist = []
    undonelist = []
    wronglist = []
    
    with open(samplelist, 'r') as indata:
        for line in indata:
            lline = line.strip().split('\t')
            lib = lline[0]
            cell = lline[1]
            name = lline[3]
            path = lline[2]
            
            querypath = os.path.join('/', path, 'final_summary.txt')
           
            if os.path.exists(querypath):
                donelist.append(name)
                cmd = "mkdir -p {projdir}/RawData/{name}/{lib}/{cell}/ &&\\\n".format(**locals())
                cmd += "cat {path}/fastq_pass/*fastq | pigz > {projdir}/RawData/{name}/{lib}/{cell}/{cell}.fastq.gz &&\\\n".format(**locals())
                cmd += "cp {path}/*sequencing_summary.txt {projdir}/RawData/{name}/{lib}/{cell}/{cell}.sequencing_summary.txt &&\\\n".format(**locals())
                
                with open((name + '.' +suffix + '.sh'), 'w') as outdata:
                    outdata.write(cmd)
                    
            elif not os.path.exists(querypath) and os.path.exists(os.path.join('/', path)):
                undonelist.append(name)
                   
            else:
                wronglist.append([name, path])
                
    return (donelist, undonelist, wronglist)            
    

            
if __name__ == "__main__":
    samplelist = sys.argv[1]
    projdir = sys.argv[2]
    donelist, undonelist, wronglist = writeshell(samplelist, projdir) 
    
    ##print done sample:
    
    for name in donelist:
        print(add_color('\t'.join((str(done), ">>>")), colorama.Fore.GREEN)),
        print(add_color(' %s     拆分完成，开始生成合并数据脚本.' % name, colorama.Fore.RED) )
        done += 1
    print('\n\n')
    ## print undone sample:
    for name in undonelist:
        print(add_color('\t'.join((str(undone), ">>>")), colorama.Fore.GREEN)),
        print(add_color(' %s 未拆分完成， be patient.' % name, colorama.Fore.YELLOW) )
        undone += 1
        
    for info in wronglist:
        print(add_color('\n%s\tthe path is invalid, please check !!\n \t%s\t%s' % (wrong, info[0], info[1]), colorama.Fore.BLUE))
        wrong += 1
    
            