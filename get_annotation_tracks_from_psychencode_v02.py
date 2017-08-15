#!/usr/bin/python

##########


import re
import glob
import os
import sys, argparse

prog_name = 'get_annotation_tracks_from_psychencode_v02.p'

#anno="/home/akhan/projects/non_coding_variants/iMEGES/layer1/gwava_unmatched_TN.avinput"

def main():
    parser = argparse.ArgumentParser(description='The intergenic genes', prog = prog_name)
    parser.add_argument('-f1', '--anno',  required = True, metavar = 'annovar file', type = str, help ='The annovar file')


    args = parser.parse_args()
    anno = os.path.abspath(args.anno)


    dir="H3K27Ac"
    os.chdir(dir)

    def get_annotate(anno):
#    dir="H3K27Ac"
 #   os.chdir(dir)
        for file in glob.glob("*.bed"):
            file11= file.rsplit( ".", 1 )[ 0 ]  
            os.system('''annotate_variation.pl '''+ anno +''' ~/projects/non_coding_variants/feature_vector_annotation/humandb -bedfile '''+file+''' --buildver hg19 -dbtype bed -regionanno -out '''+file11+'''"_bed" ''')

#if __name__ == '__main__':
    get_annotate(anno)

    for file in glob.glob("*.hg19_bed"):
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

    os.chdir("/home/akhan/projects/non_coding_variants/iMEGES/layer1")
#######

    dir="H3K27me3"
    os.chdir(dir)
    def get_annotate(anno):
  #  dir="H3K27me3"
   # os.chdir(dir)
        for file in glob.glob("*.bed"):
            file11= file.rsplit( ".", 1 )[ 0 ]
            os.system('''annotate_variation.pl '''+ anno +''' ~/projects/non_coding_variants/feature_vector_annotation/humandb -bedfile '''+file+''' --buildver hg19 -dbtype bed -regionanno -out '''+file11+'''"_bed" ''')

	get_annotate(anno)

	
    for file in glob.glob("*.hg19_bed"):
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
    os.chdir("/home/akhan/projects/non_coding_variants/iMEGES/layer1")
################
    
    
    dir="H3K4me1"
    os.chdir(dir)
    def get_annotate(anno):
#    dir="H3K4me1"
 #   os.chdir(dir)
        for file in glob.glob("*.bed"):
            file11= file.rsplit( ".", 1 )[ 0 ]
            os.system('''annotate_variation.pl '''+ anno +''' ~/projects/non_coding_variants/feature_vector_annotation/humandb -bedfile '''+file+''' --buildver hg19 -dbtype bed -regionanno -out '''+file11+'''"_bed" ''')

#if __name__ == '__main__':
    get_annotate(anno)

    for file in glob.glob("*.hg19_bed"):
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
    os.chdir("/home/akhan/projects/non_coding_variants/iMEGES/layer1")

#################

    dir="H3K4me3"
    os.chdir(dir)
    def get_annotate(anno):
#    dir="H3K4me3"
 #   os.chdir(dir)
        for file in glob.glob("*.bed"):
            file11= file.rsplit( ".", 1 )[ 0 ]
            os.system('''annotate_variation.pl '''+ anno +''' ~/projects/non_coding_variants/feature_vector_annotation/humandb -bedfile '''+file+''' --buildver hg19 -dbtype bed -regionanno -out '''+file11+'''"_bed" ''')

#if __name__ == '__main__':
    get_annotate(anno)

    for file in glob.glob("*.hg19_bed"):
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

if __name__ == '__main__':
    main()
