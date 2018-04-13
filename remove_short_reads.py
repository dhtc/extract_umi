#!/usr/bin/env python
__author__='dh'
__version__='v1.0'
import sys
def remove_short_reads(fq):
    h_tol=[]
    seq_tol=[]
    h2_tol=[]
    sc_tol=[]
    with open(fq,'r')as f:
        h=f.readline()
        seq=f.readline()
        h2=f.readline()
        sc=f.readline()
        while h:
            if len(seq)>30:
                h_tol.append(h)
                seq_tol.append(seq)
                h2_tol.append(h2)
                sc_tol.append(sc)
            h=f.readline()
            seq=f.readline()
            h2=f.readline()
            sc=f.readline()
    with open(fq[:-2]+'rm_sh.fq','w')as f:
        for n in range(len(h_tol)):
            f.write(h_tol[n])
            f.write(seq_tol[n])
            f.write(h2_tol[n])
            f.write(sc_tol[n])
if __name__=='__main__':
    fq=sys.argv[1]
    remove_short_reads(fq)
