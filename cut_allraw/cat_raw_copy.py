#!/usr/bin/python
#-*- coding:utf-8 -*-

"""
120例线下数据，由于一个样本的存在多条lane的情况，先合并，在分析；
"""

import sys
import glob
import os

sam_list = []
with open('name', 'r') as indata:
    for line in indata:
        line = line.strip()
        sam_list.append(line)
        
        

def cat_raw(sam, raw_path):
    import glob
    sam = sam
    sam_path = os.path.join(raw_path, sam)
    print(sam_path)
    read1_path = os.path.join(raw_path, sam, '*_1.fq.gz')
    reads1_list = glob.glob(read1_path)
    #reads2_list = glob.glob(read2_path)      ##在某些情况下会报错
    reads2_list = [x.replace('_1.fq.gz', '_2.fq.gz') for x in reads1_list]
    cmd1 = 'zcat {reads1_list} | gzip -c > {sam_path}/{sam}_L1_1.fq.gz '.format(reads1_list = ' '.join(reads1_list), sam_path = sam_path, sam = sam )
    cmd2 = 'zcat {reads2_list} | gzip -c > {sam_path}/{sam}_L1_2.fq.gz '.format(reads2_list = ' '.join(reads2_list), sam_path = sam_path, sam = sam )
    
    with open(os.path.join(sam_path, 'merge_raw1.sh'), 'w') as inf1:
        inf1.write(cmd1)
    with open(os.path.join(sam_path, 'merge_raw2.sh'), 'w') as inf2:
        inf2.write(cmd2)
    
    
def main():

    for sam in sam_list:
        cat_raw(sam, raw_path)
        shell = 'sh {0} &&\\'.format(raw_path + '/' + sam + '/merge_raw1.sh') 
        shell += '\n sh {0} '.format(raw_path + '/' + sam + '/merge_raw2.sh') 
    
        with open(os.path.join(raw_path, sam, 'cat_all.sh'), 'w') as out:
            out.write(shell)
        cmd = 'qsub -V -l vf=4g  {0}'.format(os.path.join(raw_path, sam, 'cat_all.sh'))
        os.system(cmd)
        

if __name__ == '__main__':
    
    raw_path = sys.argv[1]
    main()


