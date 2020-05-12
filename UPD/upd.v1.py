#!/usr/bin/python
#coding: utf-8
#Last Modified : Thu Jun 25 16:33:48 CST 2015
#Modify : output line's gene name form AAchange
#Modify : recognise format .:.:.:.
#Modify 20160818 : add check family phenotype, at least one patient per family
#Modify 20160908 : add check depth of multi-allele mutation to make sure compound heterozygous
#Modify 20171011 : add check compound heterozygous, remove site in all normal samples.
#Modify 20170424 : debug X dominant for male when 0/0 in file
#Modify 20180119: recessive het het -> hom, 患者父母必须是杂合或者没覆盖，其他同胞是非纯合
#Modify 20180528: 之前增加过一次纯合突变的患者，其父母必须要是杂合的，不能是0/0，但是没有考虑男性X染色体，纯合的时候，其父亲可以是0/0；本次修正这个问题
#Modify 20180614: 同一行，有多个基因满足复合杂合时，不重复输出行
#Modify 20191218: 当孩子为P时，才将父母的名称添加至parelation列，此外，有家系关系的样本，必须添加上父母列
import os
import sys
import re
import argparse

from argparse import RawTextHelpFormatter
parser = argparse.ArgumentParser(description = 'Function: Dominant, Recessive and Compound heterozygous model analysis.Notice that all gene name in AAChange column must in GeneName column\nContact: yuhuan@novogene.com', formatter_class = RawTextHelpFormatter)
#parser.add_argument('-M', '--model', metavar = 'String', required = True, help = 'Dominant(D)|Recessive(R)|Compound heterozygous(C).')
parser.add_argument('-S', '--samp_info', metavar = 'File', required = True, help = '#FamilyID, SampleID and Normal/Patient must be inclued, samples in the same family must have the same familyID.')
#parser.add_argument('-Fid', '--FamilyID',metavar = 'String', required = True, help = 'in sampleinfo')
parser.add_argument('-I', '--in', metavar = 'File', required = True, help = 'The merged file use to do model analysis with .xls(.gz) format.')
parser.add_argument('-T', '--type',metavar = 'String', required = True, help = 'Mutation type, snp,indel,or both',choices=['snp','indel','snp.indel'])
parser.add_argument('-AC', '--allcase', metavar = 'String', required = False, help = 'Mutations are same in all case, use for sporadic sample, Y for all case must have same mutaion, N for case can be ref. choices=[Y, N], defalut=Y (yes)', choices=['Y', 'N'], default='Y')
parser.add_argument('-O', '--out', metavar = 'file', default  = './', help = 'The output files prefix')
argv = vars(parser.parse_args())

#model     = argv['model'].strip()
samp_info = argv['samp_info'].strip()
infile    = argv['in'].strip()
type      = argv['type'].strip()
out       = argv['out'].strip()
#Fid       = argv['FamilyID'].strip()
allcase   = argv['allcase'].strip()
#define the open file handle
def safe_open(file_name,mode):
        try:
                if not file_name.endswith('.gz'):
                        return open(file_name,mode)
                else:
                        import gzip
                        return gzip.open(file_name,mode)
        except IOError:
                print file_name + ' do not exist!'

def get_het(sample):
    if re.match(r'\d+',sample):
        sample_format = re.search('(\d+)[\/|\|](\d+)',sample)
        if sample_format.group(1) == sample_format.group(2) and sample_format.group(1) != '0':
            return 'hom'
        elif sample_format.group(1) == sample_format.group(2) and sample_format.group(1) == '0':
            return 'ref'
        elif sample_format.group(1) != sample_format.group(2) and \
            (sample_format.group(1) == '0' or sample_format.group(2) == '0'):
            return 'het'
        else:
            return 'chet'  ##compound heterozygous
    else:
        return 'nocall'  

def getpos2(name,title,name2=''):
    ntitle = [i.lower() for i in title]
    if name.lower() in ntitle:
        pos = ntitle.index(name.lower())
    elif name2 != '' and name2.lower() in ntitle:
        pos = ntitle.index(name2.lower())
    else:
        if name2 != '':
            exit('Warning: %s and %s not in title.' %(name,name2))
        else:
            exit('Warning: %s not in title.' %name)
    return pos

def getRatio(formatinfo,dpinfo):
    dpdic = dict(zip(formatinfo,dpinfo))
    if ('DV' in dpdic) and re.match(r'\d+',dpdic['DV']) and re.match(r'\d+',dpdic['DP']) and int(dpdic['DP']) > 0:
        dpratio = float(dpdic['DV'])/int(dpdic['DP'])
    elif ('AD' in dpdic) and  re.match(r'\d+',dpdic['AD']) and re.match(r'\d+',dpdic['DP']) and int(dpdic['DP']) > 0:
        Adp = 0
        for adp in range(1,len(dpdic['AD'].split(','))):
            if dpdic['AD'].split(',')[adp]!=".":
                Adp += float(dpdic['AD'].split(',')[adp])
        dpratio = float(Adp)/int(dpdic['DP'])
    else:
        dpratio = 1
        print "Please check the FORMAT %s %s cloum, if AD/DV and DP in this colum or if DP value larger than 0, dpratio set as 1." %(formatinfo, dpinfo)
    return dpratio,dpdic


def get_upd(gt_proband,gt_father,gt_mother):#gt_proband,gt_father.gt_mother
    tag = []#UN/BI/UA/UI pa/ma
    if get_het(gt_proband) in ['ref','hom']:
        if get_het(gt_father) in ['ref','hom'] and get_het(gt_mother) in ['ref','hom'] and gt_father != gt_mother:
            tag.append('UA')
            if gt_proband == gt_father:
                tag.append('pa')
            elif gt_proband == gt_mother:
                tag.append('ma')
        elif (get_het(gt_father) in ['ref','hom']) and (gt_proband != gt_father) and (get_het(gt_mother) == 'het'):
            tag.append('UI')
            tag.append('ma')
        elif (get_het(gt_mother) in ['ref','hom']) and (gt_proband != gt_mother) and (get_het(gt_father) == 'het'):
            tag.append('UI')
            tag.append('pa')
        else:
            tag.append('UN')
    elif get_het(gt_proband) == 'het':
        if get_het(gt_father) in ['ref','hom'] and get_het(gt_mother) in ['ref','hom'] and gt_father != gt_mother:
            tag.append('BI')
        else:
            tag.append('UN')
    elif (get_het(gt_proband) in ['ref','nocall']) and (get_het(gt_father) in ['ref','nocall']) and (get_het(gt_mother) in ['ref','nocall']):
        tag.append('noVariant')
    else:
        tag.append('UN')
    #print '_'.join(tag)
    return '_'.join(tag)


def main():
    with open (samp_info,'r') as f1, safe_open(infile,'r') as f2, open(out+'.xls','w') as f3, open(out+".stat.xls",'w') as f4:
        famD = {}
        statD = {}#statistic for the UN ,UA, UI, BI, percent
        #get the sample's relation
        print ">>>>>>>>>>Start to get the information ...."
        for i in f1:
            ii = i.strip()
            if ii.startswith("#") and 'familyid' not in ii.lower() :
                continue
            elif ii.lower().startswith('#familyid') and ('pa' in ii.lower()) and ('ma' in ii.lower()):
                ii_h = ii.lower().split('\t')#the header
                #print ii_h
            else:
                iiline = ii.split('\t')
                #print iiline
                if len(iiline) == len(ii_h):#the line contine the proband and the pa,ma
                    famid = iiline[ii_h.index('#familyid')]
                    probandid = iiline[ii_h.index('sampleid')]
                    paid = iiline[ii_h.index('pa')]
                    maid = iiline[ii_h.index('ma')]
                    if not famD.get(famid,False):
                        famD.setdefault(famid,[]).append(probandid)
                        famD.setdefault(famid,[]).append(paid)
                        famD.setdefault(famid,[]).append(maid)
                        #statD[famid] = {'total':0,'UN':0,'UA':0,'UI':0,'BI':0}
                        statD[famid] = {}
                        for chr in [i for i in range(1,22+1)]:
                            statD[famid][str(chr)] = {'total':0,'UN':0,'UA':0,'UI':0,'BI':0}
                    else:
                        print "The family has not the only one proband!Please check! "+ famid + "\t" + probandid
                else:
                    print "The proband %s has not the information of his father and mother. " %iiline[ii_h.index('sampleid')]
        print famD
        print statD 
        print ">>>>>>>>>>>>>>>Start judge:........."
        #judge the variants for infile
        for line in f2:
            iline = line.strip().split('\t')
            #print iline[0]
            if iline[0].startswith("##"):
                continue
            elif iline[0] in ("Priority","#CHROM"):
                headerlist = iline
                print headerlist
                pformat = line.strip().lower().split("\t").index('format')
                pChr = getpos2('#chrom',line.strip().lower().split("\t"),'chrom')
        #        print pChr
                outheader = ['tag','UPDfam'] + headerlist
                f3.write("\t".join(outheader)+"\n")
                for famid in famD:
                    for iisam in famD[famid][0:4]:
#                        famD.setdefault(famid,[]).append(headerlist.index(iisam))##add the index for the sample to get the genotype
                        famD[famid].append(headerlist.index(iisam))
                print famD
                
            else:
                gtIndex = iline[pformat].split(":").index("GT")
                tag = [] #UN,UA,UI,BI
                updfam = [] #only the UA and UI familyID
                chrom = iline[pChr]
                #print chrom
                if chrom not in ('X','Y'):
                    for famid in famD:
                        gt_proband = iline[famD[famid][3]].split(":")[gtIndex]
                        gt_father = iline[famD[famid][4]].split(":")[gtIndex]
                        gt_mother = iline[famD[famid][5]].split(":")[gtIndex]
                        upd = get_upd(gt_proband,gt_father,gt_mother)#UA_pa/ma,UI_pa/ma,BI,UN,noVariant
                 #   print upd+": "+gt_proband+" "+gt_father+" "+gt_mother
                        tag.append(upd+"_"+famid)
                        statD[famid][chrom]['total'] += 1
                        statD[famid][chrom][upd.split("_")[0]] += 1
                        if upd not in ("UN","noVariant","BI"):
                            updfam.append(famid)
                    if updfam != []:
                        outlist = [';'.join(tag),';'.join(updfam)] + iline
                    else:
                        outlist = [';'.join(tag),'--'] + iline
                    f3.write("\t".join(outlist)+"\n")
        print statD
        f4.write("\t".join(["FamID","CHR","Total","UI","UIperc","UA","UAperc","BI","BIperc","UN","UNperc"])+"\n")
        for famid in famD:
            for ichr in [str(i) for i in range(1,22+1)]:
                if statD[famid][ichr]['total'] != 0:
                    uiPerc = statD[famid][ichr]['UI']/float(statD[famid][ichr]['total'])
                    uaPerc = statD[famid][ichr]['UA']/float(statD[famid][ichr]['total'])
                    biPerc = statD[famid][ichr]['BI']/float(statD[famid][ichr]['total'])
                    unPerc = statD[famid][ichr]['UN']/float(statD[famid][ichr]['total'])
                    statout = '%s\t%s\t%d\t%d\t%.4f\t%d\t%.4f\t%d\t%.4f\t%d\t%.4f' %(famid,ichr, statD[famid][ichr]['total'], statD[famid][ichr]['UI'], uiPerc, statD[famid][ichr]['UA'], uaPerc,statD[famid][ichr]['BI'],biPerc,statD[famid][ichr]['UN'],unPerc)
#            f4.write("\t".join(statout)+"\n")
                else:
                    statout = '%s\t%s\t0\t0\t-\t0\t-\t0\t-' %(famid,ichr)
                f4.write(statout+"\n")
if __name__ == '__main__':
    main()                    
          



