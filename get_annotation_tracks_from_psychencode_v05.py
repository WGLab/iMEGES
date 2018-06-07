#!/usr/bin/env python

##########


from __future__ import division
import re
import glob
import os
import sys, argparse


###################


    
def get_annotate1(bed1, annPath):
    
    print ("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! NoTICE: Running for H3K27Ac !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ")
    for file in glob.glob("H3K27Ac/*.bed"):
        file1=file.rsplit("/")[1]
        file2= file.rsplit( ".", 1 )[0]  
        #os.system('''annotate_variation.pl '''+ bed1 +''' humandb -bedfile '''+file1+''' --buildver hg19 -dbtype bed -regionanno -out '''+file2+'''"_bed" ''')
        os.system(annPath+'''/annotate_variation.pl '''+ bed1 +''' H3K27Ac -bedfile '''+file1+''' --buildver hg19 -dbtype bed -regionanno -out '''+file2+'''"_bed" ''')

    for file in glob.glob("H3K27Ac/*.hg19_bed"):
        f = open(( file.rsplit( ".", 1 )[ 0 ] ) + ".txt", "w")
        file1=open(file,'r')
        for j in file1.readlines():
            j = j.strip()
            j1=j.split('\t')
            j2=j1[2]
            if (j2=="chrX") or (j2=="chrY"):
                j3=j2
                hg19='b37'
                j4=str(j3)+"_"+str(j1[3])+"_"+str(j1[5])+"_"+str(j1[6])+"_"+str(hg19)
                f.write('%s\n' %j4)
            else:
                j5=j2
                hg19='b37'
                j6=str(j5)+"_"+str(j1[3])+"_"+str(j1[5])+"_"+str(j1[6])+"_"+str(hg19)
                f.write('%s\n' %j6)
    f.close()
    print ("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Notice: Done for H3K27Ac !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    
    #os.system("rm *_bed_bed *_bed.hg19_bed *_bed.log")
if __name__ == '__main__':
    get_annotate1(bed1)



#######

def get_annotate2(bed1, annPath):
    print ("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! NOTICE: running for H3K27me3 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    for file in glob.glob("H3K27me3/*.bed"):
        file1=file.rsplit("/")[1]
        file2= file.rsplit( ".", 1 )[ 0 ]
        #os.system('''annotate_variation.pl '''+ bed1 +''' humandb/ -bedfile '''+file1+''' --buildver hg19 -dbtype bed -regionanno -out '''+file2+'''"_bed" ''')
        os.system(annPath+'''/annotate_variation.pl '''+ bed1 +''' H3K27me3 -bedfile '''+file1+''' --buildver hg19 -dbtype bed -regionanno -out '''+file2+'''"_bed" ''')


	
    for file in glob.glob("H3K27me3/*.hg19_bed"):
        f = open(( file.rsplit( ".", 1 )[ 0 ] ) + ".txt", "w")
        file1=open(file,'r')
        for j in file1.readlines():
            j = j.strip()
            j1=j.split('\t')
            j2=j1[2]
            if (j2=="chrX") or (j2=="chrY"):
                j3=j2
                hg19='b37'
                j4=str(j3)+"_"+str(j1[3])+"_"+str(j1[5])+"_"+str(j1[6])+"_"+str(hg19)
                f.write('%s\n' %j4)
            else:
                j5=j2
                hg19='b37'
                j6=str(j5)+"_"+str(j1[3])+"_"+str(j1[5])+"_"+str(j1[6])+"_"+str(hg19)
                f.write('%s\n' %j6)
    f.close()
    print ("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! NOTICE: Done for H3K27me3 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
if __name__ == '__main__':
    get_annotate2(bed1)



################

    
def get_annotate3(bed1, annPath):
    print ("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! NOTICE: Running for H3K4me1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    for file in glob.glob("H3K4me1/*.bed"):
        file1=file.rsplit("/")[1]
        file2= file.rsplit( ".", 1 )[ 0 ]
        #os.system('''annotate_variation.pl '''+ bed1 +''' humandb/ -bedfile '''+file1+''' --buildver hg19 -dbtype bed -regionanno -out '''+file2+'''"_bed" ''')
        os.system(annPath+'''/annotate_variation.pl '''+ bed1 +''' H3K4me1 -bedfile '''+file1+''' --buildver hg19 -dbtype bed -regionanno -out '''+file2+'''"_bed" ''')


    for file in glob.glob("H3K4me1/*.hg19_bed"):
        f = open(( file.rsplit( ".", 1 )[ 0 ] ) + ".txt", "w")
        file1=open(file,'r')
        for j in file1.readlines():
            j = j.strip()
            j1=j.split('\t')
            j2=j1[2]
            if (j2=="chrX") or (j2=="chrY"):
                j3=j2
                hg19='b37'
                j4=str(j3)+"_"+str(j1[3])+"_"+str(j1[5])+"_"+str(j1[6])+"_"+str(hg19)
                f.write('%s\n' %j4)
            else:
                j5=j2
                hg19='b37'
                j6=str(j5)+"_"+str(j1[3])+"_"+str(j1[5])+"_"+str(j1[6])+"_"+str(hg19)
                f.write('%s\n' %j6)
    f.close()
    print ("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! NOTICE: Done for H3K4me1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")

if __name__ == '__main__':
    get_annotate3(bed1)

#################

    
def get_annotate4(bed1, annPath):
    print ("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! NOTICE: Running for 3K4me3 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    for file in glob.glob("H3K4me3/*.bed"):
        file1=file.rsplit("/")[1]
        file2= file.rsplit( ".", 1 )[ 0 ]
        #os.system('''annotate_variation.pl '''+ bed1 +''' humandb/ -bedfile '''+file1+''' --buildver hg19 -dbtype bed -regionanno -out '''+file2+'''"_bed" ''')
        os.system(annPath+'''/annotate_variation.pl '''+ bed1 +''' H3K4me3 -bedfile '''+file1+''' --buildver hg19 -dbtype bed -regionanno -out '''+file2+'''"_bed" ''')


    for file in glob.glob("H3K4me3/*.hg19_bed"):
        f = open(( file.rsplit( ".", 1 )[ 0 ] ) + ".txt", "w")
        file1=open(file,'r')
        for j in file1.readlines():
            j = j.strip()
            j1=j.split('\t')
            j2=j1[2]
            if (j2=="chrX") or (j2=="chrY"):
                j3=j2
                hg19='b37'
                j4=str(j3)+"_"+str(j1[3])+"_"+str(j1[5])+"_"+str(j1[6])+"_"+str(hg19)
                f.write('%s\n' %j4)
            else:
                j5=j2
                hg19='b37'
                j6=str(j5)+"_"+str(j1[3])+"_"+str(j1[5])+"_"+str(j1[6])+"_"+str(hg19)
                f.write('%s\n' %j6)
    f.close()
    print ("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! NOTICE: Done for 3K4me3 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
if __name__ == '__main__':
    get_annotate4(bed1)






