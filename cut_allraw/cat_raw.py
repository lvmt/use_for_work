#!/usr/bin/python
#-*- coding:utf-8 -*-

"""
适用于广符儿1100例，数据释放
"""

import sys
import glob
import os


def get_pos(name, title):
    ntitle = [x.lower() for x in title ]
    if name.lower() in ntitle:
        pos = ntitle.index(name.lower())
    else:
        print('name is not in title')
        exit()
    return pos
    
def get_sam(sam_info):
    sam_list  = []
    with open(sam_info, 'r') as indata:
        for line in indata:
            line = line.strip()
            lline = line.split('\t')
            if line.startswith('#') and len(lline)  < 5:
                pass 
            elif line.startswith('#') and len(lline) >= 5 :
                sampos = get_pos('sampleid', lline)
            else:
                sam = lline[sampos]
                sam_list.append(sam)
    return sam_list

def cat_raw(sam, raw_path):
    import glob
    sam = sam
    sam_path = os.path.join(raw_path, sam)
    print(sam_path)
    read1_path = os.path.join(raw_path, sam, '*_1.fq.gz')
    #read2_path = os.path.join(raw_path, sam, '*_2.fq.gz')   ##在某些情况下会报错
    reads1_list = glob.glob(read1_path)
    #reads2_list = glob.glob(read2_path)      ##在某些情况下会报错
    reads2_list = [x.replace('_1.fq.gz', '_2.fq.gz') for x in reads1_list]
    cmd1 = 'zcat {reads1_list} | gzip -c > {sam_path}/{sam}_read1.raw.fq.gz '.format(reads1_list = ' '.join(reads1_list), sam_path = sam_path, sam = sam )
    cmd2 = 'zcat {reads2_list} | gzip -c > {sam_path}/{sam}_read2.raw.fq.gz '.format(reads2_list = ' '.join(reads2_list), sam_path = sam_path, sam = sam )
    md5 = '''
    cd {sam_path}  &&\\
        md5sum {sam}_read1.raw.fq.gz > {sam}_read1.raw.fq.gz.MD5.txt  && \\
        md5sum {sam}_read2.raw.fq.gz > {sam}_read2.raw.fq.gz.MD5.txt'''.format(**locals())

    with open(os.path.join(sam_path, 'merge_raw1.sh'), 'w') as inf1:
        inf1.write(cmd1)
    with open(os.path.join(sam_path, 'merge_raw2.sh'), 'w') as inf2:
        inf2.write(cmd2)
    with open(os.path.join(sam_path, 'merge_md5.sh'), 'w') as inf3:
        inf3.write(md5)

    
def main():
    sam_list = get_sam(sam_info)
    for sam in sam_list:
        cat_raw(sam, raw_path)
        shell = 'sh {0} &&\\'.format(raw_path + '/' + sam + '/merge_raw1.sh') 
        shell += '\n sh {0} &&\\'.format(raw_path + '/' + sam + '/merge_raw2.sh') 
        shell += '\n sh {0}'.format(raw_path + '/' + sam + '/merge_md5.sh')
        with open(os.path.join(raw_path, sam, 'cat_all.sh'), 'w') as out:
            out.write(shell)
        cmd = 'qsub -V -l vf=4g  {0}'.format(os.path.join(raw_path, sam, 'cat_all.sh'))
        os.system(cmd)
        

if __name__ == '__main__':
    sam_info = sys.argv[1]
    raw_path = sys.argv[2]
    main()


