#!/usr/bin/env python
from __future__ import division
__author__='dh'
__version__='v2.20'
'''
use for extract umi from sequence and add to header
'''
import os,re,gzip,sys
import numpy as np

if '2.' not in sys.version:
    print('based on python2')
'''
global function
'''

def get_file(fqs_file_path):

    with open(fqs_file_path,'r')as f:
        fqs=[n.strip('\n').split()for n in f.readlines()]
    return fqs
class extract_umi(object):
    def __init__(self,fqs,pattern='.{0,8}(.{8})T{3}'):
        self.fqs=fqs
        self.pattern=pattern
    # @staticmethod
    # def word_count(strings,list):
    #     wordlist=list
    #     wordcount={}
    #     for string in strings:
    #         if word not in wordcount:
    #             wordcount.setdefault(word,1)
    #         else:
    #             wordcount[word]+=1
    #     return wordcount
    def extract_umi(self):
        pattern=re.compile(self.pattern)
        umi={}
        umi_uniq={}
        umi_uniq2=[]
        for fq in self.fqs:
            fq1=fq[0]
            fq2=fq[1]
            h_tol=[]
            seq_tol=[]
            h2_tol=[]
            sc_tol=[]
            fq1_h = []
            fq1_seq = []
            fq1_h2 = []
            fq1_sc = []
            if fq2.endswith('gz'):
                pass
            elif fq2.endswith(('fq','fastq')):
                    while fq1_h:
                        fq1_h.append(f.readline())
                        fq1_seq.append(f.readline())
                        fq1_h2.append(f.readline())
                        fq1_sc.append(f.readline())
                with open(fq2,'r')as f:
                    with open(fq1, 'r')as f1:
                        h=f.readline()
                        seq=f.readline()
                        h2=f.readline()
                        sc=f.readline()
                        h1=f1.readline()
                        seq1=f1.readline()
                        h21=f1.readline()
                        sc1=f1.readline()
                        while h:
                            res=pattern.search(seq)
                            if res:
                                h_tol.append(' '.join([h.split(' ')[0]+'_'+res.group(1),h.split(' ')[1]]))
                                seq_tol.append(seq[res.end():].lstrip('T'))
                                h2_tol.append('+\n')
                                sc_tol.append(sc[-len(seq[res.end():].lstrip('T')):])
                                fq1_h.append(' '.join([h1.split(' ')[0]+'_'+res.group(1),h1.split(' ')[1]]))
                                fq1_seq.append(seq1)
                                fq1_h2.append(h21)
                                fq1_sc.append(sc1)
                                umi.setdefault(fq,[])
                                umi[fq].append(res.group(1)+'\n')
                                # finished=(fqs.index(fq)+1)/len(fqs)*()
                                # process_bar='['+'#'*finished+' '*(1-finished)+']'+'%.2f'%finished+'%'+'\r'
                                # sys.stdout.write(process_bar)
                                # sys.stdout.flush()
                            h=f.readline()
                            seq=f.readline()
                            h2=f.readline()
                            sc=f.readline()
                with open(fq2[:-3]+'_extract.fq','w')as f:
                    for n in range(len(h_tol)):
                        f.write(h_tol[n])
                        f.write(seq_tol[n])
                        f.write(h2_tol[n])
                        f.write(sc_tol[n])
                    print '%s finished'%fq2
                with open(fq1[:-3]+'_extract.fq','w')as f:
                    for n in range(len(fq1_h)):
                        f.write(fq1_h[n])
                        f.write(fq1_seq[n])
                        f.write(fq1_h2[n])
                        f.write(fq1_sc[n])
                    print '%s finished'%fq1
                umi_np=np.array(umi[fq2])
                umi_uniq.setdefault(fq2,np.unique(umi_np))
                print "%s count"%fq2
                t=len(umi_uniq[fq2])
                count=1
                print "start count umi in %s"%fq2
                for n in umi_uniq[fq2]:
                    umi_uniq2.append("%s  %d\n"%(n.strip('\n'),len(np.where(umi_np==n)[0])))
                    finished=count/t
                    count+=1
                    process_bar='['+'#'*int(finished*80)+' '*int((1-finished)*80)+']'+'%.2f'%(finished*100)+'%'+'\r'
                    sys.stdout.write(process_bar)
                    sys.stdout.flush()
                with open(fq2[:-3]+'_umi_uniq.txt','w')as f:
                    f.writelines(umi_uniq2)
            else:
                print 'must be fasq(fq) or fq.gz'
        with open('umi_sum.txt','w')as f:
            umi_sum_out=[]
            umi_sum=np.array([m for n in umi_uniq for m in n])
            umi_uniq_sum=np.unique(umi_sum)
            umi_uniq_count_sum=[len(np.where(umi_sum==n)[0])for n in umi_uniq_sum]
            for u,n in zip(umi_uniq_sum,umi_uniq_count_sum):
                umi_sum_out.append("%s  %d\n"%(u,n))
            f.writelines(umi_sum_out)


if __name__=='__main__':
    fqs_file_path=sys.argv[1]
    fqs=get_file(fqs_file_path)
    #pattern=sys.argv[2]
    ext=extract_umi(fqs)
    ext.extract_umi()
