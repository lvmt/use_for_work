#/usr/bin/env python
#-*- coding:utf-8 -*-

"""检测ONT下机数据
对多个flowcell或者同个flowcell运行多次的情况进行判断
"""

import os
import sys
import time
from collections import defaultdict
sys.path.append('/ifs/TJPROJ3/DISEASE/share/Software/python/site-packages')
import colorama


colorama.init()
suffix = time.strftime("%Y_%m_%d", time.localtime())


## get pos
def getpos(name, title):
    ntitle = [x.lower() for x in title]
    if name.lower() in ntitle:
        pos = ntitle.index(name.lower())
    else:
        print("%s is not in %s" % (name.lower(), ntitle))
        
    return pos

## 给拆分状态添加颜色标识
def add_color(text, color=colorama.Fore.BLUE):
    return '%s%s%s' % (color, text, colorama.Fore.RESET)

## 设置背景颜色
def backcolor(text, color=colorama.Back.CYAN):
    return '%s%s%s' % (color, text, colorama.Back.RESET)

## 得到样本对应文库，flowcell及其下机路径字典信息
def getdict(infile):
    alldict = defaultdict(dict)
    """
    {'sam1':{'ont':
                {'cell1':['path1'], 'cell2':['path1', path2]}
            }
        
    }
    """
    with open(infile, 'r') as indata:
        for line in indata:
            line = line.strip()
            if line.startswith('#'):
                title = line.strip('#').split('\t')
                lib = getpos('libid', title)
                cell = getpos('flowcell', title)
                path = getpos('path', title)
                sam = getpos('sample', title)
            else:
                line = line.strip().split('\t')
                s_lib = line[lib]
                s_cell = line[cell]
                s_path = line[path]
                s_sam = line[sam]
                
                if not s_sam in alldict.keys():
                    alldict[s_sam][s_lib] = defaultdict(dict)
                    alldict[s_sam][s_lib][s_cell] = defaultdict(set)
                    alldict[s_sam][s_lib][s_cell] = {s_path}
                
                else:
                    if s_cell in alldict[s_sam][s_lib].keys():
                        alldict[s_sam][s_lib][s_cell].add(s_path)
                    else :
                        alldict[s_sam][s_lib][s_cell] = defaultdict(set)
                        alldict[s_sam][s_lib][s_cell] = {s_path}
                        
    return  alldict


## 判断路径是否存在，或者是否拆分完成
def judge(query):
    Flag = ""
    queryfile = os.path.join(query, 'final_summary.txt')  ## summary 存在，表明拆分完成
    if os.path.exists(query):
        if os.path.exists(queryfile):
            Flag = "done"
        else:
            Flag = "doning"
    else:
        Flag = "bad path"
        
    return Flag

doneName = []      ## 记录完全拆分好的样本名称
undoneName = []    ## 记录没有完全拆分好的样本名称
##  获取所有样本的详细信息
def getdetail(infile):
    """
    先判断是否拆分完成，如果拆分完成，先合并
    """
    alldict = getdict(infile)
    allsam = alldict.keys()
    
    n = 1
    for sam in allsam:
        print(backcolor('%s、' % n, colorama.Back.YELLOW) ), 
        print(sam)
        ont = alldict[sam].keys()[0]
        print('  >>>%s' % ont)
        flowcell = alldict[sam][ont]    ##  返回flowcell + path 的字典
        """
        defaultdict(<type 'dict'>, {'PAE02015': set(['/TJPROJ5/DISEASE/ONT_rawdata/FZTD19H001532-1A/OND00525/20191004_0811_2-A9-D9_PAE02015_54b06977']), 
        'PAE05474': set(['/TJPROJ5/DISEASE/ONT_rawdata/FZTD19H001532-1A/OND00525/20190930_1100_2-A5-D5_PAE05474_fe909c2e'])})
        """
        samtmpdict = defaultdict(dict)  ## 记录整个样本是否完全拆分完成
        """{cell1:status, cell2:status}
        """
        for cell,path in flowcell.items():       # 一个样本可能存在多个cell
            samtmpdict[cell] = defaultdict(list)
            tmpdict = defaultdict(list)         ## 记录每个样本的每个flowcell 的拆分状态
            print('    >>>%s\t%s'% (cell, path))
            for p in path:                      # 一个cell可能存在多条路径  
                    Flag = judge(p)
                    tmpdict[Flag].append(p)
                    
            if len(tmpdict) == 1 and  "done" in tmpdict:
                print(add_color('\t****%s 拆分完成 \n' % cell, colorama.Fore.RED))
                samtmpdict[cell] = "done"
            else:
                for f,p in tmpdict.items():
                    print(add_color('\t****%s,%s\n' % (f, p), colorama.Fore.BLUE) )
                samtmpdict[cell] = "undone"
                
        sam_tag = {x for x in samtmpdict.values()}
        if len(sam_tag) == 1 and "done" in sam_tag:
            print(backcolor('\tall done !!!!!!!!\n', colorama.Back.RED))
            doneName.append(sam)
        else:
            print(backcolor('\tnot all done !!!!!!!!\n', colorama.Back.LIGHTBLACK_EX))
            undoneName.append(sam)
        
        n += 1
    return (doneName,undoneName)
        

## 按照样本名称，生成合并代码,输出合并脚本，命名：Sam + 日期 + sh
def writecmd(sam, projdir):
    cmd = ""
    alldict = getdict(infile)    
    tmpdict = alldict[sam]  
    ont = tmpdict.keys()[0]
    flowcell = tmpdict[ont]    # {cell:[path], cell2:[path, path3]}  组成的字典

    ontdir = os.path.join(projdir, 'RawData', sam,  ont)
    os.system('mkdir -p %s' % os.path.join(projdir, 'RawData'))
              
    with open(os.path.join(projdir, 'RawData', sam + '.' +suffix + '.sh'), 'w') as shell:
        for cell,path in flowcell.items():
            celldir = os.path.join(ontdir, cell)
            cmd = "############ %s\nmkdir -p %s  &&\\\n\n" % (cell, celldir)
            #print(cell)
            ## 合并fastq文件
            cmd += "cat "
            for p in path:
                cmd += os.path.join(p, 'fastq_pass', '*fastq') + " "
            cmd += " | pigz > %s &&\\\n\n" % os.path.join(celldir, cell + '.fastq.gz')
            
            ## 合并summary文件
            cmd += "cat "
            for p in path:
                cmd += os.path.join(p, '*sequencing_summary.txt') + " "
            cmd += "> %s \n" % os.path.join(celldir, cell + '.sequencing_summary.txt')
                
            #print(cmd)
            shell.write(cmd)
        
    
def main(infile, projdir):
    #alldict = getdict(infile)
    doneName,undoneName = getdetail(infile)
    print(add_color('输入样本数：%s' % str(len(doneName) + len(undoneName) ), colorama.Fore.RED))
    print(add_color('开始生成合并原始数据脚本', colorama.Fore.GREEN))
    print(add_color('拆分完成的样本数：%s' % len(doneName), colorama.Fore.LIGHTMAGENTA_EX))
    n = 1
    for sam in doneName:
        print(add_color('%s、开始生成样本：%s 的合并脚本 ' % (n, sam), colorama.Fore.LIGHTRED_EX))
        writecmd(sam, projdir)
        print(add_color('   >>>样本：%s 合并脚本已经生成' % sam, colorama.Fore.LIGHTBLUE_EX))
        n += 1
        
        
    
    


if __name__ == "__main__":
    infile = sys.argv[1]
    projdir = '/TJPROJ5/DISEASE/Proj/X101SC19081135-Z01_shandaEye_42li_ONT/TGS.X101SC19081135.shandaEye_42li_B1'
   
    main(infile, projdir)
    
