#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   demo.py
@Time    :   2021/11/15 16:11:56
@Author  :   lvmt 
@Version :   1.0
'''

# here put the import lib
import yaml



yaml_info = yaml.load(open('vep_function.yaml'))
print(yaml_info['Priority'])