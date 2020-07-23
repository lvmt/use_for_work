#!/bin/usr/python
import os
import argparse
import time


class RemoveDup(object):

        def __init__(self,args):
            self.args = args
            self.infile = args['infile']
            self.outfile = args['outfile']
            self.dupfile = args['id']

        def safe_open(self,filename,mode='r'):
            if not filename.endswith('.gz'):
                return open(filename,mode)
            elif filename.endswith('.gz'):
                import gzip
                return gzip.open(filename,mode)
            else:
                exit('something is wrong,check..{filename}'.format(**locals()))

        def start(self):
            id_dic = {}
            with self.safe_open(self.dupfile) as fr:
                for i in fr:
                        id='@'+i.strip()
                        id_dic[id] = i

            with self.safe_open(self.infile, 'r') as fr,\
                self.safe_open(self.outfile, 'w') as fw:
                try:
                    while True:
                        a1=fr.next().strip()
                        a2=fr.next().strip()
                        a3=fr.next().strip()
                        a4=fr.next().strip()
                        fq_id=a1.split()[0]
                        if id_dic.get(fq_id):
                            #print id_dic[fq_id]
                            continue
                        else:
                            fw.write(a1+'\n'+a2+'\n'+a3+'\n'+a4+'\n')
                except StopIteration:
                    pass
 

 
def main():
    
    dd = RemoveDup(args)
    dd.start()


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='去除指定ID的dup序列')
    parser.add_argument('--id', help='Duplication ID list',required=True)
    parser.add_argument('--infile', help='输入fq.gz', required=True)
    parser.add_argument('--outfile', help='输出fq,未压缩格式', required=True)

    args = vars(parser.parse_args())
    main()
