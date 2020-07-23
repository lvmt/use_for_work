#!/usr/bin/env python
#-*- coding:utf-8 -*-


import os 
import pysam 
#import re 
#import socket
import textwrap


'''
修改去除dup脚本中的问题
'''

class RemoveDup(object):

    def __init__(self,args):
        self.args = args
        self.samplelist = args['samplelist']
        self.projdir = args['projdir']
        self.freq = args['freq']
        self.outdir = args['outdir']
        
        self.version = args['version']
        self.queues = args['queues'] 

        

    def mkdir(self,dirname):
        if not os.path.exists(dirname):
            os.system('mkdir -p {dirname}'.format(**locals()))

    def check_dir(self,dirname):
        if os.path.exists(dirname):
            return True
        return False

    def safe_open(self,filename,mode='r'):
        if not filename.endswith('.gz'):
            return open(filename,mode)
        elif filename.endswith('.gz'):
            import gzip
            return gzip.open(filename,mode)
        else:
            exit('somt thing is wrong, check!')

    def write_shell(self,cmd,shellname):
        with open(shellname,'w') as fw:
            fw.write(cmd)

    def get_dup_lane(self):
        '''处理samplelist文件,得到每行的信息
        '''
        removedup_list = []  # 记录单个样本的全部的dupid信息

        # 索引信息
        lane_index = 0
        sam_index = 2
        lib_index = 3
        path_index = 6 if self.version == 'Y' else 5

        # 处理samplelist文件
        with self.safe_open(self.samplelist, 'r') as fr:
            for line in fr:
                linelist = line.strip().split('\t')
                if line.startswith('#'):
                    headlist = [item.lower() for item in linelist]
                    if 'ori_lane' in headlist:
                        lane_index = headlist.index('ori_lane')
                    else:
                        lane_index = headlist.index('lane')
                    sam_index = headlist.index('samid')
                    path_index = headlist.index('path')
                    lib_index = headlist.index('libid')
                    continue
                sam = linelist[sam_index]
                path = linelist[path_index]
                lane = 'L' + linelist[lane_index][-1]  # L3
                lib = linelist[lib_index]

                flowcell = path.strip('/').split('/')[-1]
                lane_num = lane[-1]  # 3

                line_info = [sam,lib,lane,path,lane_num,flowcell]
                if line_info not in removedup_list:
                    removedup_list.append(line_info)
        
        return removedup_list

    def start(self):
        '''记录样本dup信息
        从finalbam中记录的，所有的lane均合并在一起啦
        所有此次提取的dupid是该样本全部文库全部lane的合集        
        '''
        removedup_list = self.get_dup_lane()
        bam_list = []

        for line_info in removedup_list:
            
            sam,lib,lane,path,lane_num,flowcell = line_info

            fq1 = '{path}/{lib}/{lib}_{lane}_1.fq.gz'.format(**locals())
            fq2 = '{path}/{lib}/{lib}_{lane}_2.fq.gz'.format(**locals())
            adapter1 = '{path}/{lib}/{lib}_{lane}_1.adapter.list.gz'.format(**locals())
            adapter2 = '{path}/{lib}/{lib}_{lane}_2.adapter.list.gz'.format(**locals())
            out_fq1 = '{self.outdir}/RemoveDup/{flowcell}/{lib}/{lib}_{lane}_1.fq'.format(**locals())
            out_fq2 = '{self.outdir}/RemoveDup/{flowcell}/{lib}/{lib}_{lane}_2.fq'.format(**locals())

            bam = '{self.projdir}/Mapping/{sam}.{sam}/{sam}.nodup.bam'.format(**locals())
            if not self.check_dir(bam):
                exit('Bam file not exists. {sam}'.format(**locals()))
                
            sam_outdir = '{self.outdir}/RemoveDup/{flowcell}/{lib}'.format(**locals())
            self.mkdir(sam_outdir)

            # 记录每个样本的dupid总信息
            if bam not in bam_list:
                bam_list.append(bam)
                dup_file = '{sam_outdir}/{sam}.DupListOri.xls'.format(**locals()) #
                bamfile = pysam.Samfile(bam, 'rb')
                with self.safe_open(dup_file, 'w') as fw:
                    for read in bamfile.fetch():
                        # skip the read is not primary. get dup only from primary reads
                        if  read.flag & 0x900 != 0:
                            continue
                        if read.is_duplicate:
                            dup_id = read.qname
                            fw.write(dup_id + '\n')
                bamfile.close()

            # 得到每个样本每条lane的dup_id文件
            udup_list = '{sam_outdir}/{sam}.DupLi_{lane}.xls'.format(**locals())
            cmd = """
            sed -n '/:{lane_num}:/'p {dup_file} > {udup_list}
            """.format(**locals())
            os.system(cmd)

            # 得到指定百分比的dup数据
            tmp_lane_dup_list = []
            with self.safe_open(udup_list,'r') as fr:
                for line in fr:
                    line = line.strip()
                    tmp_lane_dup_list.append(line)

            sort_lane_dup_list = sorted(list(set(tmp_lane_dup_list)))
            length_dup = int(len(sort_lane_dup_list))

            # 得到需要的dup数量
            use_length = int(length_dup * self.freq)
            use_dup_list = sort_lane_dup_list[0:use_length] 
            freq_list = '{sam_outdir}/{sam}.FiDupLi_{lane}.xls'.format(**locals())

            with open(freq_list,'w') as fw:
                fw.write('\n'.join(use_dup_list))

            ## write shell
            # F_shell = '{sam_outdir}/{sam}_{flowcell}_{lib}_{lane}_ReDup_1.sh'.format(**locals())
            # R_shell = '{sam_outdir}/{sam}_{flowcell}_{lib}_{lane}_ReDup_2.sh'.format(**locals())
            F_shell = '{sam_outdir}/{sam}_{lane}_ReDup_1.sh'.format(**locals()) # 因此从目了进行了区分,
            R_shell = '{sam_outdir}/{sam}_{lane}_ReDup_2.sh'.format(**locals())

            cmdF = textwrap.dedent('''
            python /ifs/TJPROJ3/DISEASE/share/Disease/RmDup/rm_dup.py \\
                --id {freq_list} \\
                --infile {fq1} \\
                --outfile {out_fq1} &&\\

            pigz -p4 -f {out_fq1}
            ''').format(**locals())

            cmdR = textwrap.dedent('''
            python /ifs/TJPROJ3/DISEASE/share/Disease/RmDup/rm_dup.py \\
                --id {freq_list} \\
                --infile {fq2} \\
                --outfile {out_fq2} &&\\

            pigz -p4  -f {out_fq2}
            ''').format(**locals())

            self.write_shell(cmdF,F_shell)
            self.write_shell(cmdR,R_shell)

            ## run shell
            cmd = textwrap.dedent('''
            # remove dup
            qsub -V -l vf=5g -q disease2.q
            ''')

def main():

    dup = RemoveDup(args)
    dup.start()


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='\033[1;32muse for removedup || 解决同lane问题\033[0m')
    parser.add_argument('--samplelist', help='记录样本的sample list文件',required=True)
    parser.add_argument('--projdir', help='项目路径',required=True)
    parser.add_argument('--version', help='项目版本',default='Y')
    parser.add_argument('--freq', help='去除的频率',default='0.5',type=float)
    parser.add_argument('--queues', help='分析队列信息')
    parser.add_argument('--outdir', help='输出目录',required=True)

    args = vars(parser.parse_args())
    main()


