#!/usr/bin/python
# -*- coding: utf-8 -*-
import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from email.header import Header
from email.utils import parseaddr, formataddr
import textwrap
from collections import defaultdict
from collections import OrderedDict
import os
import re
import jinja2
from jinja2 import PackageLoader,Environment

env = Environment(loader=PackageLoader('result', 'templates'))
template = env.get_template('demo3.html') 


class Check(object):

    def __init__(self, qcfile, mapfile, teacher):
        self.qcfile = qcfile
        self.mapfile = mapfile
        self.qc_dict = self.get_qc_dict()
        self.map_dict = self.get_mapping_dict()
        self.teacher = teacher # 林云婷 和 李娜 对clean的要求不一样


    def qc(self):
        """
        #qc
        clean: 林云婷：12G，李娜10G
        error: 小于0.1
        q30: 大于90
        q20: 大于95
        """
        clean = {}
        error = {}
        q30 = {}
        q20 = {}

        qc_failed_dict ={}
        for sam in self.qc_dict:
            if self.teacher == "lin":
                if self.qc_dict[sam]['clean'] < 12:
                    clean[sam] = self.qc_dict[sam]['clean']
            elif self.teacher == "na":
                if self.qc_dict[sam]['clean'] < 10:
                    clean[sam] = self.qc_dict[sam]['clean']

            if self.qc_dict[sam]['error'] > 0.1:
                error[sam] = self.qc_dict[sam]['error']
            if self.qc_dict[sam]['q30'] < 90:
                q30[sam] = self.qc_dict[sam]['q30']
            if self.qc_dict[sam]['q20'] < 95:
                q20[sam] = self.qc_dict[sam]['q20']
        
        qc_failed_dict['clean'] = clean
        qc_failed_dict['error'] = error
        qc_failed_dict['q30'] = q30
        qc_failed_dict['q20'] = q20
        # print(qc_failed_dict)
        return qc_failed_dict


    def mapping(self):
        """
        #map
        dup: 小于30
        20x: 大于98
        map: 大于99 mapped
        effective: 55, 捕获效率
        depth: 100x
        """
        dup = {}
        _20x = {}
        _map = {}
        effective = {}
        depth = {}
        map_failed_dict = {}
        for sam in self.map_dict:
            if float(self.map_dict[sam]['dup']) > 30:
                dup[sam] = self.map_dict[sam]['dup']
            if float(self.map_dict[sam]['20x']) < 98:
                _20x[sam] = self.map_dict[sam]['20x']
            if float(self.map_dict[sam]['mapped']) < 99:
                _map[sam] = self.map_dict[sam]['mapped']
            if float(self.map_dict[sam]['effective']) < 55:
                effective[sam] = self.map_dict[sam]['effective']
            if float(self.map_dict[sam]['depth']) < 100:
                depth[sam] = self.map_dict[sam]['depth']
        
        map_failed_dict['dup'] = dup
        map_failed_dict['20x'] = _20x
        map_failed_dict['map'] = _map
        map_failed_dict['effective'] = effective
        map_failed_dict['depth'] = depth 

        return map_failed_dict

    
    def get_qc_dict(self):
        qc_dict = defaultdict(dict)
        with open(self.qcfile, 'r') as f:
            title_list = f.readline().strip().split('\t')
            for line in f:
                linelist = line.strip().split('\t')
                sam = linelist[0]
                if sam in qc_dict:
                    qc_dict[sam]['rawdata'].append(linelist[4])
                    qc_dict[sam]['effective'].append(linelist[6])
                    qc_dict[sam]['error'].append(linelist[7])
                    qc_dict[sam]['q20'].append(linelist[8])
                    qc_dict[sam]['q30'].append(linelist[9])
                    qc_dict[sam]['gc'].append(linelist[10])
                else:
                    qc_dict[sam]['rawdata'] = [linelist[4]]
                    qc_dict[sam]['effective'] = [linelist[6]]
                    qc_dict[sam]['error'] = [linelist[7]]
                    qc_dict[sam]['q20'] = [linelist[8]]
                    qc_dict[sam]['q30'] = [linelist[9]]
                    qc_dict[sam]['gc'] = [linelist[10]]
        # merge all lane of one sample
        qc_merge_dict = defaultdict(dict)
        for sam in qc_dict:
            rawdata_sum = 0
            effective_sum = 0
            error_sum = 0
            q20_sum = 0
            q30_sum = 0
            gc_sum = 0
            for i in range(len(qc_dict[sam]['gc'])):
                rawdata_sum += float(qc_dict[sam]['rawdata'][i])
                effective_sum += float(qc_dict[sam]['effective'][i]) * float(qc_dict[sam]['rawdata'][i])
                error_sum += float(qc_dict[sam]['error'][i]) * float(qc_dict[sam]['rawdata'][i])
                q20_sum += float(qc_dict[sam]['q20'][i]) * float(qc_dict[sam]['rawdata'][i])
                q30_sum += float(qc_dict[sam]['q30'][i]) * float(qc_dict[sam]['rawdata'][i])
                gc_sum += float(qc_dict[sam]['gc'][i]) * float(qc_dict[sam]['rawdata'][i])
            
            merge_rawdata = rawdata_sum
            merge_effective = round(effective_sum / merge_rawdata, 2)
            merge_error = round(error_sum / merge_rawdata, 2)
            merge_q20 = round(q20_sum / merge_rawdata, 2)
            merge_q30 = round(q30_sum / merge_rawdata, 2)
            merge_gc = round(gc_sum / merge_rawdata, 2)

            qc_merge_dict[sam]['rawdata'] = merge_rawdata
            qc_merge_dict[sam]['effective'] = merge_effective
            qc_merge_dict[sam]['clean'] = round(merge_effective * merge_rawdata / 100, 2)
            qc_merge_dict[sam]['error'] = merge_error
            qc_merge_dict[sam]['q20'] = merge_q20
            qc_merge_dict[sam]['q30'] = merge_q30
            qc_merge_dict[sam]['gc'] = merge_gc

        # print(qc_merge_dict)
        return qc_merge_dict
    

    @staticmethod
    def clean_str(str):
        new_list = [x for x in re.split('\(|\)', str) if x]
        return new_list[-1].strip("%")
    

    def get_mapping_dict(self):
        mapping_dict = defaultdict(dict)
        with open(self.mapfile, 'r') as f:
            title = f.next()
            try:
                while f:
                    # print(f.next())
                    linelist = f.next().strip().split("\t")
                    mapping_dict[linelist[1]]['dup'] = Check.clean_str(linelist[3])
                    mapping_dict[linelist[1]]['mapped'] = Check.clean_str(linelist[4])
                    mapping_dict[linelist[1]]['proper_map'] = Check.clean_str(linelist[5])
                    mapping_dict[linelist[1]]['pe_map'] = Check.clean_str(linelist[6])
                    mapping_dict[linelist[1]]['effective'] = Check.clean_str(linelist[12])
                    mapping_dict[linelist[1]]['depth'] = Check.clean_str(linelist[14])
                    mapping_dict[linelist[1]]['coverage'] = Check.clean_str(linelist[16])
                    mapping_dict[linelist[1]]['20x'] = Check.clean_str(linelist[19])
                    mapping_dict[linelist[1]]['10x'] = Check.clean_str(linelist[20])
                    mapping_dict[linelist[1]]['4x'] = Check.clean_str(linelist[21])
            except StopIteration as e:
                print(e)
        # print(mapping_dict)
        return mapping_dict



# def sendMail(body, sub, to_mail):
#     """
#     body: 邮件主题
#     sub：邮件主题
#     to_mail：邮件接收人, '11.com,22.com'
#     """
#     smtp_server = 'smtp.163.com'
#     from_mail = '13554221497@163.com'
#     mail_pass = 'IEVEYITYZYDPDEDC'

#     msg = MIMEText(body, "html", "utf-8")
#     sub = "【加测】".decode("utf-8") + sub
#     msg['From'] = from_mail
#     msg['To'] = to_mail
#     msg['Subject'] = Header(sub, 'utf-8').encode()

#     try:
#         s = smtplib.SMTP()     
#         s.connect(smtp_server, "25")   
#         s.login(from_mail, mail_pass)
#         s.sendmail(from_mail, to_mail, msg.as_string())  # as_string()把MIMEText对象变成str     
#         s.quit()
#     except smtplib.SMTPException as e:
#         print("Error: %s" % e)


def get_teacher(pnfile):
    teacher_map = {
        "P101SC17101594-01": "lin",
        "P101SC17040027-01": "lin",
        "X102SC19051571-Z01": "na",
        "X101SC19080128-Z01": "na"
    }
    
    with open(pnfile, 'r') as fr:
        linelist = fr.readline().strip().split("\t")
        fenqi_num = linelist[0]
        projname = linelist[1].decode("utf-8")
        projnum = "-".join(fenqi_num.split("-")[0:2])

    teacher = teacher_map[projnum]
    return fenqi_num,projname, projnum, teacher



if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(description="check广妇儿项目问题")
    parser.add_argument("-p", "--projdir", help="项目路径")
    parser.add_argument("-j", "--jobname", help="job名称")

    args = vars(parser.parse_args())
    print(args)

    ## other params
    filepath = os.path.join(args['projdir'], 'Report', args['jobname'], 'Mapping')
    qcfile = filepath + '/qcstat.xls'
    mapfile = filepath + '/mappingstat.xls' 
    pnfile = os.path.join(args['projdir'], 'pn.txt')
    
    fenqi_num,projname, projnum, teacher = get_teacher(pnfile)

    ## check 质控 mapping 
    cc = Check(qcfile, mapfile, teacher)


    # 开始渲染html文档
    
    # qc check内容
    qc_failed_dict = cc.qc()
    cleandata = qc_failed_dict['clean']
    q30 = qc_failed_dict['q30']
    q20 = qc_failed_dict['q20']
    error = qc_failed_dict['error']

    # map check 内容
    map_failed_dict = cc.mapping()
    dup = map_failed_dict['dup']
    effective = map_failed_dict['effective']
    _map = map_failed_dict['map']
    _20x = map_failed_dict['20x']
    depth = map_failed_dict['depth'] 


    body = template.render(projname_str=projname,
                           projdir=args["projdir"],
                           fenqi_str=fenqi_num,
                           cleandata_dict=cleandata,
                           depth_dict=depth,
                           q30_dict = q30,
                           q20_dict = q20,
                           error_dict = error,
                           dup_dict=dup,
                           _20x_dict=_20x,
                           jobname=args["jobname"],
                           projnum=projnum,
                           teacher=teacher)

    

    outfile = os.path.join(args["projdir"],'Report', args['jobname'] ,'auto_check.html') 
    with open(outfile, 'w') as fw:
        fw.write(body.encode("utf-8"))

    subject = "【加测详情】"+projname.encode("utf-8")

    cmd = '''
    source ~/.bash_profile
    sendemail "lvmengting4480@novogene.com,wangshan@novogene.com,zhangruiwen@novogene.com"  "{subject}" "{outfile}" html
    '''.format(**locals())
    
    
    os.system(cmd)
    #sendMail(body, projname,'lvmengting4480@novogene.com, wangshan@novogene.com')

    
