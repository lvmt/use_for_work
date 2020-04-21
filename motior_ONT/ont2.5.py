#/usr/bin/env python
#-*- coding:utf-8 -*-

"""检测ONT下机数据
对多个flowcell或者同个flowcell运行多次的情况进行判断
"""


"""增加控制选项：
1、输出合并完成的样本名称
2、输出没有合并完成的样本名称
3、 是否生成合并脚本
4、查询每个样本的路径文库信息
add:20191012
1、生成指定样本的路径时， check指定样本是否已经拆分完成，
add:20191014
同一个flowcell合并时，需要去除summary.txt中的重复表头信息，不然质控会报错
add:20191107
生产存在数据拆分问题，在合并数据前，需要check summary文件中的行数 和 fastq_pass中的行数是否一致
"""

import os
import sys
import time
from collections import defaultdict
sys.path.append('/ifs/TJPROJ3/DISEASE/share/Software/python/site-packages')
import colorama
import argparse

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
        for cell,path in flowcell.items():                        # 一个样本可能存在多个cell
            samtmpdict[cell] = defaultdict(list)
            tmpdict = defaultdict(list)                           ## 记录每个样本的每个flowcell 的拆分状态
            print('    >>>%s\t%s'% (cell, path))
            
            for p in path:                                             # 一个cell可能存在多条路径  
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
def writecmd(sam, projdir,alldict):
    cmd = ""
    alldict = alldict  
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
            cmd += " | sed '2,${/filename_fastq/d}'> %s \n" % os.path.join(celldir, cell + '.sequencing_summary.txt')
                
            #print(cmd)
            shell.write(cmd)
        
def query(sam, alldict):
    samdict = alldict[sam]
    ont = samdict.keys()[0]
    flowcell = samdict[ont]
    print(add_color('查询样本名：%s '% sam, colorama.Fore.GREEN))
    print(add_color('  >>>ONT编号：%s '% ont, colorama.Fore.LIGHTGREEN_EX))
    
    samtmpdict = defaultdict(dict)  ## 记录整个样本是否完全拆分完成
    """{cell1:status, cell2:status}
    """
    for cell,path in flowcell.items():                        # 一个样本可能存在多个cell
        samtmpdict[cell] = defaultdict(list)
        tmpdict = defaultdict(list)                           ## 记录每个样本的每个flowcell 的拆分状态
        print('    >>>%s\t%s'% (cell, path))
        
        for p in path:                                             # 一个cell可能存在多条路径  
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
    else:
        print(backcolor('\tnot all done !!!!!!!!\n', colorama.Back.LIGHTBLACK_EX))
    

## 在输出指定样本的合并脚本时， 需要check指定的样本是否已经拆分完成
def checksam(samlist, donelist):
    wronglist = []
    for sam in samlist:
        if sam in donelist:
            pass
        else:
            wronglist.append(sam)
    
    return wronglist
   

## 针对准备分析的样本，check下summary文件中的行数及名称和fastq_pass中的是否一致
def check_summary(inpath):
    
    pass 



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="用于三代数据的处理")
    parser.add_argument('--infile', help="samplelist, 包含4列信息，libid,flowcellid,sam,path,需要#号开头", required=True)
    parser.add_argument('--detail', help="展示所有样本的拆分状态", action='store_true')
    parser.add_argument('--no', help="在打印信息的时候，是否生成脚本，", action="store_true")
    parser.add_argument('--cmd', help="生成合并脚本，默认对所有拆分完成的样本，需要先指定detail参数")   ## 默认对所有拆分完成的生成脚本，可以自己指定 
    parser.add_argument('--query', help="查询输入样本的下机信息")
    
    args = parser.parse_args()
    print(args)
    
    ##  get var
    infile = args.infile
    #projdir = args.projdir
    sam = args.query

    alldict = getdict(infile)
    
    ## 查看所有样本的拆分信息，打印简单的统计信息
    if args.detail:
        doneName,undoneName = getdetail(infile)
        print(add_color('\n样本信息统计\n输入样本数：%s' % str(len(doneName) + len(undoneName) ), colorama.Fore.RED))
        print(add_color('拆分完成的样本数：%s' % len(doneName), colorama.Fore.LIGHTMAGENTA_EX))
        for sam in doneName:
            print(sam),
        print('\n')
        print(add_color('未拆分完成的样本数：%s' % len(undoneName), colorama.Fore.LIGHTGREEN_EX))
        for sam in undoneName:
            print(sam),
            
        print(add_color('\n你可以使用query， 查看具体样本的拆分情况\n', colorama.Fore.LIGHTRED_EX))
        
        
        ##  指定cmd参数，不指定的话，默认对已经拆分完全的样本生成合并脚本
        ##  指定cmd参数，只对指定的样本，生成合并脚本
        
        if not args.no:
            print(backcolor('%s选择生成合并脚本%s' % ('*' * 6, '*' * 6), colorama.Back.YELLOW))
            projdir = input("\n请输入项目路径：")
            if args.cmd:
                samlist = args.cmd.split(',')
                print(add_color('生成**指定样本**的合并脚本：', colorama.Fore.YELLOW))
                print(add_color('  >>指定样本名称：%s' % ' '.join(samlist), colorama.Fore.LIGHTGREEN_EX))
                
                ## check 样本信息
                wronglist = checksam(samlist, doneName)
                if len(wronglist) == 0:
                    print(add_color('  >>>指定样本全都拆分完成', colorama.Fore.RED))
                else:
                    print(add_color('    >>>please check!!! the sample 【 %s 】 还未拆分完成，请query！！' % " ".join(wronglist), colorama.Fore.LIGHTGREEN_EX))
                    exit(add_color('please check first!!!', colorama.Fore.LIGHTRED_EX))
                
                ## write shell
                n = 1
                for sam in samlist:
                    print(add_color('%s、开始生成样本：%s 的合并脚本 ' % (n, sam), colorama.Fore.LIGHTRED_EX))
                    writecmd(sam, projdir, alldict)
                    print(add_color('   >>>样本：%s 合并脚本已经生成' % sam, colorama.Fore.LIGHTBLUE_EX))
                    n += 1
                
            else:
                samlist = doneName
                n = 1
                print(backcolor('生成拆分完成样本的合并脚本：', colorama.Back.YELLOW))
                for sam in samlist:
                    print(add_color('%s、开始生成样本：%s 的合并脚本 ' % (n, sam), colorama.Fore.LIGHTRED_EX))
                    writecmd(sam, projdir, alldict)
                    print(add_color('   >>>样本：%s 合并脚本已经生成' % sam, colorama.Fore.LIGHTBLUE_EX))
                    n += 1
        else:
            pass


    if args.query:
        query(sam, alldict)
   
   