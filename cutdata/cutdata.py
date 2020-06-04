#!/usr/bin/env python
#-*- coding:utf-8 -*-

import textwrap


class CutData(object):

    def __init__(self, args):
        self.sample_list = args['sample_list']
        self.target = args['target'] if args['target'] else None
        self.num = args['num'] if args['num'] else None   # 自己指定截取的数目
        self.outdir = args['outdir']
        self.info_dict = self.get_info_dict()  

    
    def get_info_dict(self):
        """
        info_dict = {
            'sam1': [path, libid, laneid],
            'sam2': [path, libid, laneid]
        }
        """

        info_dict = {}
        with open(self.sample_list, 'r') as fr:
            for line in fr:
                
                if line.startswith('#'):
                    linelist = [item.lower() for item in line.strip().split('\t')]
                    index_dict = {item:index for index,item in enumerate(linelist) }
                    print(index_dict)
                else:
                    linelist = line.strip().split('\t')
                    info_dict[linelist[index_dict['sampleid']]] = [linelist[index_dict['path']], linelist[index_dict['libid']], \
                                                                    linelist[index_dict['#ori_lane']]   ]

        return info_dict


    def out_one_shell(self, name):
        if self.target:
            target_fq = self.info_dict[self.target][0] + '/' + self.info_dict[self.target][1] + '/' + '*L' + self.info_dict[self.target][2] +'_1.fq.gz'
       
        path = self.info_dict[name][0]
        flowcell = path.strip().split('/')[-1]
        libid = self.info_dict[name][1]
        laneid = self.info_dict[name][2]

        in_fq1 = '{path}/{libid}/{libid}_L{laneid}_1.fq.gz'.format(**locals())
        in_fq2 = '{path}/{libid}/{libid}_L{laneid}_2.fq.gz'.format(**locals())
        out_fq1 = '{self.outdir}/{flowcell}/{libid}/{libid}_L{laneid}_1.fq.gz'.format(**locals())
        out_fq2 = '{self.outdir}/{flowcell}/{libid}/{libid}_L{laneid}_2.fq.gz'.format(**locals())

        outshell = '{self.outdir}/{name}_lane{laneid}.sh'.format(**locals())
        
        if self.num:
            cmd = textwrap.dedent('''
            echo start `date`
            num={self.num}
            '''.format(**locals()))
        else:
            cmd = textwrap.dedent('''
            echo start `date`
            num=`seqkit seq -j8 -n {target_fq} | wc  -l`         
            '''.format(**locals()))

        if name != self.target:
            cmd += textwrap.dedent(''' 

            seqkit head -j8 -n ${{num}} {in_fq1}  -o  {out_fq1}  &&\\

            seqkit head -j8 -n ${{num}} {in_fq2}  -o  {out_fq2}

            echo end `date`
            '''.format(**locals())) 


            with open(outshell, 'w') as ff:
                ff.write(cmd)
        

    def start(self):

        for name in self.info_dict.keys():
            self.out_one_shell(name)



if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument('sample_list', help="需要进行数据截取的sample_list文件")
    parser.add_argument('-t','--target', help="指定一个目标样本，按照这个样本的数据量进行截取，")
    parser.add_argument('outdir', help='指定截取的数据放在那里，截取之后的数据目录和之前一毛一样')
    parser.add_argument('-n','--num', help='如果没有指定目标样本，可以自己指定read数目,他的优先级高于target')

    args = vars(parser.parse_args())

    demo = CutData(args)

    demo.start()


