#!/usr/bin/env python
# -*- coding:utf-8 -*-


# number = 0
# record_info = []   # 记录所以的信息
# #把结果放在列表里面，每次的结果用append添加
# while True:
#     anser1=input('A你认罪么？认或者不认：')
#     anser2=input('B你认罪么？认或者不认：')
#     if anser1 in ('认', '不认') and anser2 in ('认', '不认'):
#         number += 1
#         record_info.append([number, anser1, anser2])
#     else:
#         print('\033[31m错误输入： 只能选择认或者不认\033[0m')
#         continue
#     if anser1 == '认' and anser2 == '认':
#         print('\033[31m各判十年\033[0m') 
#     elif anser1 == '不认' and anser2 == '认':
#         print('\033[31mA1判20年，B判1年\033[0m')
#     elif anser1 == '认' and anser2 == '不认':
#         print('\033[31mA判1年，B判20年\033[0m')
#     elif anser1 == '不认' and anser2 == '不认':
#         print('\033[31m各判三年\033[0m')
#         print('\033[31m第{}对实验者都不认罪\033[0m'.format(number))
#         break
#     else:
#         pass 
        

# #记录每一对实验者的选择
# for item in range(len(record_info)):
#     print(record_info[item])

############################################################


number = 1
record_info = []   # 记录所以的信息

while True:
    anser1=input('A你认罪么？认或者不认：')
    anser2=input('B你认罪么？认或者不认：')
    
    if anser1 == '认' and anser2 == '认':
        print('\033[31m各判十年\033[0m') 
        record_info.append([number, anser1, anser2])
        number += 1
    elif anser1 == '不认' and anser2 == '认':
        print('\033[31mA1判20年，B判1年\033[0m')
        record_info.append([number, anser1, anser2])
        number += 1
    elif anser1 == '认' and anser2 == '不认':
        print('\033[31mA判1年，B判20年\033[0m')
        record_info.append([number, anser1, anser2])
        number += 1
    elif anser1 == '不认' and anser2 == '不认':
        print('\033[31m各判三年\033[0m')
        print('\033[31m第{}对实验者都不认罪\033[0m'.format(number))
        record_info.append([number, anser1, anser2])
        number += 1
        break
    else:
        print('\033[31m错误输入： 只能选择认或者不认\033[0m')
        pass