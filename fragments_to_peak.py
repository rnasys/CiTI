#! usr/bin/env python

"""
Author: Xiaojuan Fan
Date: 2017-5-19
E-mail: fanxiaojuan@picb.ac.cn
Description: Detect reads cluster (depth >= 1)
Note: Change code when detect mature mRNA or pre-mRNA
"""

"""
Input format:
#------------------------------------------------------------------------------ 
chr1    880462    -    1
chr1    880463    -    1
chr1    880464    -    1
"""

import argparse
from numpy.lib.function_base import average
from interval import interval
import time

def createHelp():
    """
    Create the command line interface of the program.
    """
    
    epilog_string="Any bug is welcome reported to fanxiaojuan@picb.ac.cn"
    description_string='The program is going to detect reads cluster'
    parser = argparse.ArgumentParser(description=description_string,epilog=epilog_string)
    parser.add_argument('-i', '--input file', dest='fnIn', default='plasmid.mature.filter.vcf.sort', type=str,help='input vcf file')
    parser.add_argument('-i_ref', '--refgene', dest='fnIn_ref', default='refGene.txt', type=str,help='input vcf file')
    parser.add_argument('-d', '--depth', dest='d', default=1, type=int,help='minimum depth')
    parser.add_argument('-o', '--output file', dest='fnOut', default='plasmid.mature.filter.cluster.bed', type=str,help='output cluster file')
    op=parser.parse_args()
    return op

def refgene_trim(refgene_db):
    """
    Trim refgene dataset
    """
    ref_dic = {}
    for i in range(0,len(refgene_db)):
        words = refgene_db[i].strip().split('\t')
        chr_ref = words[2]
        strand_ref = words[3]
        start_trans = int(words[4])
        end_trans = int(words[5])
        try:
            ref_dic[(chr_ref,strand_ref)].append([start_trans,end_trans])
        except:
            ref_dic[(chr_ref,strand_ref)] = [[start_trans,end_trans]]
    #print ref_dic
    return ref_dic

def interval_union(interval_list):
    """
    Calculate interval union
    """
    #print interval_list
    interval_results = interval()
    for i in range(0,len(interval_list)):
        interval_results = interval(interval_list[i]) | interval_results
    return [item for item in interval_results]

def cluster_classifier(vcf_file,refgene_dic,op):
    """
    Classify cluster
    """
    cluster = []
    coverage_list = []
    block_interval = []
    pos_list = []
    pos_list_2 = []
    coverage_list_2 = []
    i = 1
    while i < len(vcf_file):
        words = vcf_file[i].strip().split('\t')
        chr_no = words[0]
        chr_no_last = vcf_file[i-1].strip().split('\t')[0]
        base_pos = int(words[1])
        base_pos_last = int(vcf_file[i-1].strip().split('\t')[1])
        strand = words[2]
        strand_last = vcf_file[i-1].strip().split('\t')[2]
        cover = int(words[3])
        #cover_last = int(vcf_file[i-1].strip().split('\t')[3])
        #print words
        KEY = (chr_no_last,strand_last)
        A = (chr_no != chr_no_last)
        B = (strand != strand_last)
        D = (base_pos - base_pos_last > 5000) #detect mature mRNA
        #D = (base_pos - base_pos_last > 1) #detect pre-mRNA
        if KEY in refgene_dic:
            for m in range(0,len(refgene_dic[KEY])):
                C = (base_pos_last >= refgene_dic[KEY][m][0] and base_pos_last <= refgene_dic[KEY][m][1] \
                and base_pos > refgene_dic[KEY][m][1])
                if C == True:
                    break
        #print C
        if base_pos == base_pos_last:
            i += 1
            continue
        else:
            if A or B or C or D: #detect mature mRNA
            #if A or B or D: #detect pre-mRNA
                #print base_pos_last,base_pos
                if pos_list != []:
                    #print pos_list
                    #print coverage_list
                    cover_max = max(coverage_list)
                    cover_min = min(coverage_list)
                    #print cover_max
                    #print base_pos_last,base_pos
                    if cover_max == 1 or cover_min >= 2:
                        #print coverage_list
                        #print base_pos_last,base_pos
                        block_interval = interval_union(pos_list)
                        block_interval = sorted(block_interval,key = lambda block_interval:block_interval[1])
                        #print block_interval
                        summit = pos_list[coverage_list.index(max(coverage_list))][1]
                        #print summit
                        cover_ave = average(coverage_list)
                        block_length = [str(int(block_interval[z][1] - block_interval[z][0])) for z in range(0,len(block_interval))]
                        block_start = [str(int(block_interval[z][0] - block_interval[0][0])) for z in range(0,len(block_interval))]
                        cluster.append('\t'.join([chr_no_last,str(int(block_interval[0][0])),str(int(block_interval[-1][-1])),'0',str(cover_ave),strand_last,\
                                                  str(int(block_interval[0][0])),str(int(block_interval[-1][-1])),'0',str(len(block_interval)),\
                                                  ','.join(block_length),','.join(block_start),str(summit),str(cover_max)]))
                    else:
                        #print coverage_list
                        #print pos_list
                        #print base_pos_last,base_pos
                        for a in range(0,len(coverage_list)):
                            if ((coverage_list[a] == 1) or (a == len(coverage_list)-1)) and pos_list_2 != []:
                                block_interval = interval_union(pos_list_2)
                                block_interval = sorted(block_interval,key = lambda block_interval:block_interval[1])
                                #print pos_list_2
                                #print base_pos_last,base_pos
                                #print block_interval
                                summit = pos_list_2[coverage_list_2.index(max(coverage_list_2))][1]
                                #print summit
                                cover_ave = average(coverage_list_2)
                                cover_max = max(coverage_list_2)
                                block_length = [str(int(block_interval[z][1] - block_interval[z][0])) for z in range(0,len(block_interval))]
                                block_start = [str(int(block_interval[z][0] - block_interval[0][0])) for z in range(0,len(block_interval))]
                                cluster.append('\t'.join([chr_no_last,str(int(block_interval[0][0])),str(int(block_interval[-1][-1])),'0',str(cover_ave),strand_last,\
                                                          str(int(block_interval[0][0])),str(int(block_interval[-1][-1])),'0',str(len(block_interval)),\
                                                          ','.join(block_length),','.join(block_start),str(summit),str(cover_max)]))
                                coverage_list_2 = []
                                block_interval = []
                                pos_list_2 = []
                            elif coverage_list[a] >= 2:
                                pos_list_2.append(pos_list[a])
                                coverage_list_2.append(coverage_list[a])
                coverage_list = []
                block_interval = []
                pos_list = []
            coverage_list.append(cover)
            pos_list.append([base_pos - 1,base_pos])#Detect mature mRNA
        i += 1
    return cluster

if __name__ == '__main__':
    time_start = time.time()
    op = createHelp()
    
    refgene_db = open(op.fnIn_ref).readlines()
    refgene_dic = refgene_trim(refgene_db)
    
    print ('Start to detect cluster...')
    vcf_file = open(op.fnIn).readlines()
    cluster_list = cluster_classifier(vcf_file,refgene_dic,op)
    print ('Cluster detection done! Items: ' + str(len(cluster_list)))
    
    out_cluster = open(op.fnOut,'w')
    out_cluster.write('\n'.join(cluster_list) + '\n')
    
    print ('Done!')
    print ('Time used: ' + str(time.time() - time_start) + ' s.')
