#! usr/bin/env python

"""
Author: Xiaojuan Fan
Date: 2016-6-5
E-mail: fanxiaojuan@picb.ac.cn
Description: Count z-score of motif
"""

"""
input file format:
#------------------------------------------------------------------------------ 
motif    count    frequency
GCGTT    292    0.0546714098483
AAATG    1374    0.257255195656
GCCCG    459    0.0859389627411
GCCCA    1019    0.190788241902
AAATC    1061    0.198651937839
"""

import argparse
import math

def createHelp():
    """
    Create the command line interface of the program.
    """
    
    epilog_string="Any bug is welcome reported to fanxiaojuan@picb.ac.cn"
    description_string='The program is going to count z-score of motif.'
    parser = argparse.ArgumentParser(description=description_string,epilog=epilog_string)
    parser.add_argument('--i-s1', '--input sample 1', dest='fnIn_s1', default='utr5.motif.frequency.txt', help='input sample 1')
    parser.add_argument('--i-s2', '--input sample 2', dest='fnIn_s2', default='utr5.bg.motif.frequency.txt', help='input sample 2 as background')
    parser.add_argument('-o', '--output file', dest='fnOut', default='utr5.motif', help='output file')
    op=parser.parse_args()
    return op

if __name__ == '__main__':
    op = createHelp()
    
    print 'Build motif dictionary in background dataset.'
    
    background = open(op.fnIn_s2).readlines()
    bg_dic = {}
    for i in xrange(0,len(background)):
        words_bg = background[i].strip().split('\t')
        motif_bg = words_bg[0]
        count_bg = int(words_bg[1])
        freq_bg = float(words_bg[2])
        bg_dic[motif_bg] = [count_bg,freq_bg]
        
    #print bg_dic
    print 'Background dic done!'
    
    print 'Start to calculate z-score...'
    final_list = []
    circRNA = open(op.fnIn_s1).readlines()
    for i in xrange(0,len(circRNA)):
        words_cir = circRNA[i].strip().split('\t')
        motif_cir = words_cir[0]
        count_cir = int(words_cir[1])
        freq_cir = float(words_cir[2])
        if motif_cir in bg_dic:
            if freq_cir != 0:
                p = (count_cir + bg_dic[motif_cir][0])/((count_cir / freq_cir + bg_dic[motif_cir][0] / bg_dic[motif_cir][1]))
                #print p
                z_score = (freq_cir - bg_dic[motif_cir][1])/math.sqrt((1/(count_cir/freq_cir) + 1/(bg_dic[motif_cir][0]/bg_dic[motif_cir][1])) * p * (1-p))
                #print z_score
                fold_change = freq_cir / bg_dic[motif_cir][1]
                final_list.append('\t'.join(words_cir)+'\t'+str(bg_dic[motif_cir][0])+'\t'+str(bg_dic[motif_cir][1])+'\t'+str(fold_change)+'\t'+str(z_score))
        else:
            #print motif_cir
            final_list.append('\t'.join(words_cir))
            
    #print final_list
    print 'Done! Start to write to output file...'
    output = open(op.fnOut,'w')
#     final_list.sort(key = lambda l:(l[1]),reverse=True)
    output.write('\n'.join(final_list)+'\n')
    
    print 'Finish!'
            
