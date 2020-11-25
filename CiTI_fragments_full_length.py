#! usr/bin/env python

"""
Author: Xiaojuan Fan
Date: 2019-9-26
Revised date: 2020-9-17
E-mail: fanxiaojuan@picb.ac.cn
Description: Extend the paired-end reads interval nucleotides and annotate the fragments
"""

"""
input bed format (*.rmdup.bed):
#------------------------------------------------------------------------------ 
chr1    100912694    100912794    2    1    +    100912694    100912794    255,0,0    1    100    0
chr1    100912751    100912832    1    1    -    100912751    100912832    255,0,0    1    81    0
chr17    79479044    79479144    1    3    +    79479044    79479144    255,0,0    1    100    0
chr17    79479310    79479381    2    3    -    79479310    79479381    255,0,0    1    71    0
chr17    39777844    39777936    1    2    +    39777844    39777936    255,0,0    1    92    0
chr17    39777915    39777976    2    2    -    39777915    39777976    255,0,0    1    61    0
"""

import argparse
import time
import portion as P

def createHelp():
    """
    Create the command line interface of the program.
    """
    
    epilog_string="Any bug is welcome reported to fanxiaojuan@picb.ac.cn"
    description_string='The program is going to extend the paired-end reads interval nucleotides'
    parser = argparse.ArgumentParser(description=description_string,epilog=epilog_string)
    parser.add_argument('-i', '--input file', dest='fnIn', default='plasmid.rmdup.bed', type=str,help='input file')
    parser.add_argument('-i_gtf', '--gtf', dest='gtf', default='gencode.v28lift37.annotation.gtf', type=str,help='input gencode gtf annotation file')
    parser.add_argument('-db', '--db', dest='db', default='gencode_v28lift37_comprehensive.txt', type=str,help='input file')
    parser.add_argument('-o', '--overhang', dest='o', default=0, type=int,help='input file')
    parser.add_argument('-o_anno', '--out-annotation', dest='anno', default='plasmid.fragment.rmdup.annotation.bed', type=str,help='output mature mRNA')
    op=parser.parse_args()
    return op

def Extract_gene_type(gtf_file):
    """
    Extract gene type from gtf file
    """
    gene_type_dic = {}
    for i in range(0,len(gtf_file)):
        if '##' not in gtf_file[i]:
            row = gtf_file[i].strip().split('\t')
            if row[2] == 'transcript':
                trans_id = row[8].split('transcript_id "')[1].split('";')[0]
                #print trans_id
                gene_type_dic[trans_id] = row[8].split('transcript_type "')[1].split('";')[0]
    return gene_type_dic

def gencode_dic(gencode_file,gene_type_dic):
    """
    Build gencode database dictionary
    """
    gen_dic = {}
    for i in range(1,len(gencode_file)):
        words_gen = gencode_file[i].strip().split('\t')
        chr_no = words_gen[2]
        trans_id = words_gen[1]
        cds_info = words_gen[13]
        cde_info = words_gen[14]
        gene_type = gene_type_dic[trans_id]
        gene_name = words_gen[12]
        TSS_start = int(words_gen[4])
        TSS_end = int(words_gen[5])
        CDS_start = int(words_gen[6])
        CDS_end = int(words_gen[7])
        strand = words_gen[3]
        start_list = [int(x) for x in words_gen[9].split(',')[:-1]]
        end_list = [int(x) for x in words_gen[10].split(',')[:-1]]
        exon_no = int(words_gen[8])
#         if (chr_no,trans_id) in gen_dic:  #Some trans_id are not unique, especially transcripts in chrX and chrY
#             print trans_id
        interval_list = [P.closedopen(start_list[x],end_list[x]) for x in range(0,exon_no)]
        interval_merge = P.empty()
        for i in range(0,len(interval_list)):
            interval_merge = interval_merge | interval_list[i]
        if gene_type == 'protein_coding':
            if (cds_info == 'cmpl') and (cde_info == 'cmpl'):
                # print (interval_merge)
                gen_dic.setdefault((chr_no,strand),[]).append([TSS_start,TSS_end,CDS_start,CDS_end,\
                                                    gene_name,gene_type,interval_merge])
        else:
            gen_dic.setdefault((chr_no,strand),[]).append([TSS_start,TSS_end,CDS_start,CDS_end,\
                                                    gene_name,gene_type,interval_merge])
    return gen_dic

def reads_collaps(bed_file):
    """
    Unique fragments
    """
    CITI_fragments = {}
    fragments_coverage = {}
    fragments_count = 0
    for i in range(0,len(bed_file),2):
        #temp_count_before = len(mature_mRNA) + len(pre_mRNA)
        left_reads = bed_file[i].strip().split('\t')
        right_reads = bed_file[i+1].strip().split('\t')
        #----------------------------------------------------- #left reads
        chr_left = left_reads[0]
        start_left = int(left_reads[1]) # start position in bed file is 0-based, but in gencode database file is 1-based
        end_left = int(left_reads[2])
        reads_cov = int(left_reads[4])
        interval_start_left = [int(x)+start_left for x in left_reads[11].split(',')]
        interval_length_left = [int(x) for x in left_reads[10].split(',')]
        interval_left = [(interval_start_left[a],interval_start_left[a]+interval_length_left[a]) for a in range(0,len(interval_start_left))]
        #print interval_left
        interval_left_reads = P.empty()
        for m in range(0,len(interval_left)):
            # print (interval_left[m])
            interval_left_reads = interval_left_reads | P.closedopen(interval_left[m][0],interval_left[m][1])
        #------------------------------------------------------# right reads
        chr_right = right_reads[0]
        start_right = int(right_reads[1])
        end_right = int(right_reads[2])
        interval_start_right = [int(x)+start_right for x in right_reads[11].split(',')]
        interval_length_right = [int(x) for x in right_reads[10].split(',')]
        interval_right = [(interval_start_right[a],interval_start_right[a]+interval_length_right[a]) for a in range(0,len(interval_start_right))]
        #print interval_right
        interval_right_reads = P.empty()
        for n in range(0,len(interval_right)):
            interval_right_reads = interval_right_reads | P.closedopen(interval_right[n][0],interval_right[n][1])
        #print interval_left_reads,interval_right_reads
#------------------------------------------------------------------------------ left reads overlapped with right reads
        if '2' in left_reads[3] and '1' in right_reads[3]:
            strand = '+'
        elif '1' in left_reads[3] and '2' in right_reads[3]:
            strand = '-'
        else:
            print ('Strand error! Please check: ' + '\t'.join(left_reads) + '\n' + '\t'.join(right_reads))
        block_start = min(start_left,start_right)
        block_end = max(end_left,end_right)
        interval_final = interval_left_reads | interval_right_reads
        if end_left > start_right:
            CITI_fragments[(chr_left,block_start,block_end,strand)] = \
                        (interval_final,'complete')
        else:
            CITI_fragments[(chr_left,block_start,block_end,strand)] = \
                        (interval_final,'incomplete')
        fragments_coverage[(chr_left,block_start,block_end,strand)] = \
                        fragments_coverage.get((chr_left,block_start,block_end,strand),0) + reads_cov
        fragments_count += reads_cov
    print ('Total fragments: ' + str(fragments_count))
    return (CITI_fragments,fragments_coverage)

def paired_interval_extend(uniq_fragment,fragment_cov,gtf_dic):
    """
    Extend the paired-end reads interval nucleotides
    """
    out_dic = {}
    total_reads = 0
    for key in uniq_fragment.keys():
        chr_no = key[0]
        #print (frag_start,frag_end)
        frag_strand = key[3]
        interval_comp = uniq_fragment[key][0]
        complete_info = uniq_fragment[key][1]
        frag_cov = fragment_cov[key]
        total_reads += frag_cov
        geneNA = 'NA'
        geneType = 'NA'
        geneRegion = 'NA'
        flag = 0
        for trans in gtf_dic[(chr_no,frag_strand)]:
            frag_start,frag_end = key[1:3]
        # for trans in gtf_dic[('chr1','-')]:
            # if chr_no == 'chr1' and frag_strand == '-':
            if frag_start > trans[0] and frag_end < trans[1]:
                #print 'Hello!'
                # print (trans)
                geneNA = trans[4]
                geneType = trans[5]
                if geneType == 'protein_coding':
                    CDS_start,CDS_end = trans[2:4]
                    if frag_start >= CDS_start and frag_end <= CDS_end:
                        geneRegion = 'CDS'
                    elif frag_strand == '+':
                        if frag_end <= CDS_start:
                            geneRegion = '5UTR'
                        elif frag_start < CDS_start and frag_end > CDS_start:
                            geneRegion = '5UTR-CDS'
                        elif frag_start < CDS_end and frag_end > CDS_end:
                            geneRegion = 'CDS-3UTR'
                        elif frag_start >= CDS_end:
                            geneRegion = '3UTR'
                    elif frag_strand == '-':
                        if frag_end <= CDS_start:
                            geneRegion = '3UTR'
                        elif frag_start < CDS_start and frag_end > CDS_start:
                            geneRegion = 'CDS-3UTR'
                        elif frag_start < CDS_end and frag_end > CDS_end:
                            geneRegion = '5UTR-CDS'
                        elif frag_start >= CDS_end:
                            geneRegion = '5UTR'
                else:
                    geneRegion = 'Null'
                # print (frag_start,frag_end,CDS_start,CDS_end,geneNA,geneRegion)
#------------------------------------------------------------------------------ intersect of fragments interval and exons interval
                frag_intersect = interval_comp & trans[-1]
                interval_comp_length = sum([interval_comp[a].upper- interval_comp[a].lower for a in range(0,len(interval_comp))])
                # print (interval_comp)
                # print (frag_intersect)
#------------------------------------------------------------------------------ fragments located in introns
                if frag_intersect == P.empty(): 
                    flag = 1
                    start_out = []
                    length_out = []
                    for interval_region in list(interval_comp):
                        start_out.append(str(int(interval_region.lower - frag_start)))
                        length_out.append(str(int(interval_region.upper - interval_region.lower)))
                    out_dic.setdefault((chr_no,frag_start,frag_end,frag_strand),[]).append((chr_no,str(frag_start),str(frag_end),\
                                            geneNA,geneType,frag_strand,\
                                            str(frag_start),str(frag_end),'intron',str(len(start_out)),\
                                            ','.join(length_out),','.join(start_out),str(frag_cov),flag,complete_info))
                else:
                    if complete_info == 'complete':
                        flag = 3
                        #print interval_comp
#------------------------------------------------------------------------------ reduce alignment noise
                        frag_intersect_length = sum([frag_intersect[a].upper-frag_intersect[a].lower for a in range(0,len(frag_intersect))])
                        absolute_diff = abs(frag_intersect_length-interval_comp_length)
                        if absolute_diff == 0:
#------------------------------------------------------------------------------ 
                            start_region = []
                            length_region = []
                            for region in frag_intersect:
                                start_region.append(str(int(region.lower - frag_start)))
                                length_region.append(str(int(region.upper - region.lower)))
                            out_dic.setdefault((chr_no,frag_start,frag_end,frag_strand),[]).append((chr_no,str(frag_start),str(frag_end),\
                                                    geneNA,geneType,frag_strand,\
                                                    str(frag_start),str(frag_end),geneRegion,str(len(start_region)),\
                                                    ','.join(length_region),','.join(start_region),str(frag_cov),flag,complete_info))
                        else:
                            start_region = []
                            length_region = []
                            for region in interval_comp:
                                start_region.append(str(int(region.lower - frag_start)))
                                length_region.append(str(int(region.upper - region.lower)))
                            out_dic.setdefault((chr_no,frag_start,frag_end,frag_strand),[]).append((chr_no,str(frag_start),str(frag_end),geneNA,geneType,\
                                                        frag_strand,str(frag_start),str(frag_end),'intron-containing',str(len(start_region)),\
                                                        ','.join(length_region),','.join(start_region),str(frag_cov),flag,complete_info))
                    else:
                        #print interval_comp
                        #print frag_intersect
#------------------------------------------------------------------------------ fragments boundaries located in exons
                        #print frag_intersect[0][0],frag_start,frag_intersect[-1][1],frag_end
                        #print abs_position
                        # print (P.closedopen(frag_start,frag_end),trans[-1])
                        interval_update = P.closedopen(frag_start,frag_end) & trans[-1]
                        # print (interval_update)
                        frag_trans_length = sum([interval_update[a].upper-interval_update[a].lower for a in range(0,len(interval_update))])
                        absolute_diff = abs(frag_trans_length-interval_comp_length)
                        #print absolute_diff
                        #print geneRegion
                        #print interval_comp
                        #print abs_position
                        if absolute_diff <= 300: #insert sequence length <=200nt
                            #print frag_trans_length,interval_comp_length
                            #print geneRegion
                            flag = 2
                            start_out = []
                            length_out = []
                            for interval_region in list(interval_update):
                                start_out.append(str(int(interval_region.lower - frag_start)))
                                length_out.append(str(int(interval_region.upper - interval_region.lower)))
                            out_dic.setdefault((chr_no,frag_start,frag_end,frag_strand),[]).append((chr_no,str(frag_start),str(frag_end),\
                                                geneNA,geneType,frag_strand,\
                                                str(frag_start),str(frag_end),geneRegion,str(len(start_out)),\
                                                ','.join(length_out),','.join(start_out),str(frag_cov),flag,complete_info))
                        else:
                            # print (trans)
                            flag = 1
                            start_out = []
                            length_out = []
                            for interval_region in list(interval_comp):
                                start_out.append(str(int(interval_region.lower - frag_start)))
                                length_out.append(str(int(interval_region.upper - interval_region.lower)))
                            out_dic.setdefault((chr_no,frag_start,frag_end,frag_strand),[]).append((chr_no,str(frag_start),str(frag_end),\
                                                geneNA,geneType,frag_strand,\
                                                str(frag_start),str(frag_end),'intron-containing',str(len(start_out)),\
                                                ','.join(length_out),','.join(start_out),str(frag_cov),flag,complete_info))
        if flag == 0:
            start_out = []
            length_out = []
            for interval_region in list(interval_comp):
                start_out.append(str(int(interval_region.lower - frag_start)))
                length_out.append(str(int(interval_region.upper - interval_region.lower)))
            out_dic[(chr_no,frag_start,frag_end,frag_strand)] = [(chr_no,str(frag_start),str(frag_end),'intergenic','intergenic',frag_strand,\
                                                 str(frag_start),str(frag_end),geneRegion,str(len(start_out)),\
                                                 ','.join(length_out),','.join(start_out),str(frag_cov),flag,complete_info)]
    print ('Total treated fragments: ' + str(total_reads))
    return out_dic
                    
def fragment_length_filter(fragment_anno_dic):
    """
    Keep the longest fragments
    """
    out_list = []
    total_fragment = 0
    for key in fragment_anno_dic.keys():
        #print fragment_anno_dic[key]
        fragments_flag = []
        fragments_length = []
        fragments_region = []
        total_fragment += int(fragment_anno_dic[key][0][-3])
        reads_coverage = [x[-3] for x in fragment_anno_dic[key]]
        if len(list(set(reads_coverage))) != 1:
            print (fragment_anno_dic[key])
        if len(fragment_anno_dic[key]) == 1:
            fragment_anno_dic[key][0] = list(fragment_anno_dic[key][0])
            fragment_anno_dic[key][0][-2] = str(fragment_anno_dic[key][0][-2])
            out_list.append('\t'.join(fragment_anno_dic[key][0]))
        else:
            for i in range(0,len(fragment_anno_dic[key])):
                fragment_anno_dic[key][i] = list(fragment_anno_dic[key][i])
                iso = fragment_anno_dic[key][i]
                iso_length = sum([int(x) for x in iso[10].split(',')])
                fragments_length.append(iso_length)
                fragments_flag.append(iso[-2])
                fragments_region.append(iso[8])
            #print fragment_anno_dic[key]
#---------------------------------------------------------------- complete fragments (Set region preference)
            region_complete = [''] * len(fragments_flag)
            max_flag = max(fragments_flag)
            #print fragments_length,fragments_region,fragments_flag
            if max_flag == 3:
                for x in range(0,len(fragments_flag)):
                    if fragments_flag[x] == max_flag:
                        fragment_anno_dic[key][x][-2] = str(fragment_anno_dic[key][x][-2])
                        region_complete[x] = fragments_region[x]
                # Set preference
                if 'CDS' in region_complete:
                    out_list.append('\t'.join(fragment_anno_dic[key][region_complete.index('CDS')]))
                elif '5UTR' in region_complete:
                    out_list.append('\t'.join(fragment_anno_dic[key][region_complete.index('5UTR')]))
                elif '3UTR' in region_complete:
                    out_list.append('\t'.join(fragment_anno_dic[key][region_complete.index('3UTR')]))
                elif '5UTR-CDS' in region_complete:
                    out_list.append('\t'.join(fragment_anno_dic[key][region_complete.index('5UTR-CDS')]))
                elif 'CDS-3UTR' in region_complete:
                    out_list.append('\t'.join(fragment_anno_dic[key][region_complete.index('CDS-3UTR')]))
                elif 'intron' in region_complete:
                    out_list.append('\t'.join(fragment_anno_dic[key][region_complete.index('intron')]))
                elif 'intron-containing'  in region_complete:
                    out_list.append('\t'.join(fragment_anno_dic[key][region_complete.index('intron-containing')]))
                elif 'Null' in region_complete:
                    out_list.append('\t'.join(fragment_anno_dic[key][region_complete.index('Null')]))
                else:
                    print (fragment_anno_dic[key])
                    print ('Gene type error!')
#----------------------------------------------------------------- incomplete fragments (choose the longest fragments)
            elif max_flag == 2:
                max_length_list = [0] * len(fragments_length)
                max_region_list = [''] * len(fragments_length)
                for y in range(0,len(fragments_flag)):
                    if fragments_flag[y] == max_flag:
                        max_length_list[y] = fragments_length[y]
                #print max_length_list
                max_length = max(max_length_list)
                #print max_length
                for z in range(0,len(max_length_list)):
                    if max_length_list[z] == max_length:
                        fragment_anno_dic[key][z][-2] = str(fragment_anno_dic[key][z][-2])
                        max_region_list[z] = fragments_region[z]
                #print max_region_list
                # Set preference
                if 'CDS' in max_region_list:
                    out_list.append('\t'.join(fragment_anno_dic[key][max_region_list.index('CDS')]))
                elif '5UTR' in max_region_list:
                    out_list.append('\t'.join(fragment_anno_dic[key][max_region_list.index('5UTR')]))
                elif '3UTR' in max_region_list:
                    out_list.append('\t'.join(fragment_anno_dic[key][max_region_list.index('3UTR')]))
                elif '5UTR-CDS' in max_region_list:
                    out_list.append('\t'.join(fragment_anno_dic[key][max_region_list.index('5UTR-CDS')]))
                elif 'CDS-3UTR' in max_region_list:
                    out_list.append('\t'.join(fragment_anno_dic[key][max_region_list.index('CDS-3UTR')]))
                elif 'intron' in max_region_list:
                    out_list.append('\t'.join(fragment_anno_dic[key][max_region_list.index('intron')]))
                elif 'intron-containing'  in region_complete:
                    out_list.append('\t'.join(fragment_anno_dic[key][region_complete.index('intron-containing')]))
                elif 'Null' in max_region_list:
                    out_list.append('\t'.join(fragment_anno_dic[key][max_region_list.index('Null')]))
            elif max_flag == 1: #Not annotated to exon region
                fragment_anno_dic[key][fragments_flag.index(1)][-2] = str(fragment_anno_dic[key][fragments_flag.index(1)][-2])
                # print (fragment_anno_dic[key])
                out_list.append('\t'.join(fragment_anno_dic[key][fragments_flag.index(1)]))
            elif max_flag == 0: #Not annotated to intragenic region
                fragment_anno_dic[key][0][-2] = str(fragment_anno_dic[key][0][-2])
                out_list.append('\t'.join(fragment_anno_dic[key][0]))
            else:
                print (fragment_anno_dic[key])
                print ('Please check flag information')
    print ('Total fragments after filtering 1: ' + str(total_fragment))
    return out_list
            
if __name__ == '__main__':
    time_start = time.time()
    op = createHelp()
    
    gtf_file = open(op.gtf).readlines()
    gene_type_dic = Extract_gene_type(gtf_file)
    
    gtf_db = open(op.db).readlines()
    gtf_dic = gencode_dic(gtf_db,gene_type_dic)
    print ('GTF dic done!')
    print ('Time used: ' + str(time.time() - time_start) + ' s.')
    #print (gtf_dic)
    
    bed_file = open(op.fnIn).readlines()
    print ('Total bed item: ' + str(len(bed_file)/2.0))
     
    uniq_fragment_dic,fragment_coverage_dic = reads_collaps(bed_file)
    # print (len(uniq_fragment_dic))
    # print (uniq_fragment_dic[('chr16',19620295,19627468,'+')])
    # print (fragment_coverage_dic[('chr16',19620295,19627468,'+')])
    print ('Time used: ' + str(time.time() - time_start) + ' s.')
     
    fragment_anno_dic = paired_interval_extend(uniq_fragment_dic,fragment_coverage_dic,gtf_dic)
    # print (fragment_anno_dic[('chr16',19620295,19627468,'+')])
    # print ('Total unique fragments: ' + str(len(fragment_anno_dic)))
    #print (fragment_anno_dic)
    print ('Time used: ' + str(time.time() - time_start) + ' s.')
     
    fragment_filter_list = fragment_length_filter(fragment_anno_dic)
    print ('Total fragments after filtering 2: ' + str(len(fragment_filter_list)))
   
    #print (fragment_filter_list)
    out_uniq_fragment = open(op.anno,'w')
    out_uniq_fragment.write('\n'.join(fragment_filter_list) + '\n')
     
    print ('Time used: ' + str(time.time() - time_start) + ' s.')
    print ('Done!')
