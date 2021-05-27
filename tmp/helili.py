#!/usr/bin/env python
#-*- coding:utf-8 -*-


'''
cnvkit  结果处理
'''

from collections import defaultdict

def safe_open(infile,mode='r'):
    if not infile.endswith('gz'):
        return open(infile, mode)
    elif infile.endswith('.gz'):
        import gzip
        return gzip.open(infile,mode)

def get_index_dict(headlist):
    index_dict = {}
    headlist = [item.lower() for item in headlist]
    for index,item in enumerate(headlist):
        index_dict[item] = index
    return index_dict

def get_info_dict(samplefile):
    '''
    info_dict 
    {
        sam:
        {'name':name,
        'bam':bam,
        'sex':sex}
    }
    '''
    info_dict = defaultdict(dict)
    with safe_open(samplefile,'r') as fr:
        for line in fr:
            linelist = line.strip().split('\t')
            sam = linelist[0]
            bam = linelist[1]
            sex = linelist[2]
            info_dict[sam] = {
                'name': sam,
                'bam': bam,
                'sex': sex
            }
    return info_dict



    pass
    

def hand_cns_file(info_dict,cnsdir):
    '''处理样本的cns文件
    '''
    for sam in info_dict:
        cns_file = '{cnsdir}/{sam}.call.cns'.format(**locals())
        sex = info_dict[sam]['sex']
        with safe_open(cns_file, 'r') as fr, open(cns_file.replace('.cns','final.cns'),'w') as fw:
            for line in fr:
                linelist = line.strip().split('\t')
                if 'start' in linelist:
                    headlist = linelist
                    fw.write('{}\n'.format('\t'.join(headlist)))
                    index_dict = get_index_dict(headlist)
                    continue
                chr = linelist[index_dict['chromosome']]
                if sex == 'female' and chr == 'Y':
                    continue
                cn = linelist[index_dict['cn']]
                if int(cn) == 2:
                    continue
                if int(cn) > 2:
                    linelist[index_dict['cn']] = 'gain'
                elif int(cn) < 2:
                    linelist[index_dict['cn']] = 'loss'
                fw.write('{}\n'.format('\t'.join(headlist)))




