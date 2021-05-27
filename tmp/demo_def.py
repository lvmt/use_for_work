#!/usr/bin/env python
#-*- coding:utf-8 -*-

import sys




infile = sys.argv[1]
outfile = sys.argv[2]




with open(infile, 'r') as fr, open(outfile + 'stat_target.xls', 'w') as fw:
    ui_and_ua = []
    all = []

    for line in fr:
        linelist = line.strip().split("_")
        if 'BI' in linelist or 'UN' in linelist:
            aa = 0
        elif 'ma' in linelist:
            aa = 1
        elif 'pa' in linelist:
            aa = -1

        if aa != 0:
            ui_and_ua.append(aa)
        
        all.append(aa)

    print(ui_and_ua)
    print(all)

    percent = len(ui_and_ua) / len(all) 

    fw.write('ui_and_ua / all : {percent}'.format(**locals()))

