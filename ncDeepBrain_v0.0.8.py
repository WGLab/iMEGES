##################################################################################
# Author: Atlas Khan (ak4046@cumc.columbia.edu)
# Created Time: 2017-05-31
# Wang genomics lab http://wglab.org/
# Description: python script for trimming of RNA-seq data using cutadapt
##################################################################################

#!/usr/bin/env python


import sys, argparse, os
from itertools import izip
import glob
import re
import numpy as np
import subprocess

from get_annotation_tracks_from_psychencode_v05 import *

from liftOver_v02 import *


version = """%prog
Copyright (C) 2017 Wang Genomic Lab
VerTect is free for non-commercial use without warranty.
Please contact the authors for commercial use.
Written by Atlas Khan, ak4046@cumc.columbia.edu and atlas.akhan@gmail.com.
============================================================================
"""

"""usage: ncDeepBrain-v06.py [-h] -f1 input.bed -E input.bed -B bPCA.R -o
                          out_name -buildver out_hg_name"""



prog_name = 'ncDeepBrain_v0.0.8.py'



def main():
    parser = argparse.ArgumentParser(description='The ncDeepBrain score for non-coding brain variants ', prog = prog_name)
    parser.add_argument('-f1', '--bed1',  required = True, metavar = 'input.bed', type = str, help ='The Annovar input/bed format file')
    parser.add_argument('-E', '--bed2',  required = True, metavar = 'input.bed', type = str, help ='The list of the known brain eqtls from CommonMind')	
    parser.add_argument('-B', '--bPCA',  required = True, metavar = 'bPCA.R', type = str, help ='The bPCA R code to do the imputation')
    parser.add_argument('-o', '--out',  required = True, metavar = 'The name of output file', type = str, help ='The name of the output file') 
    parser.add_argument('-annPath', '--table_annovar_path',  required = True, type =str, help ='The path to table_annovar and its humandb') 
    parser.add_argument('-buildver', '--hg',  required = True, metavar = 'out_hg_name', type =str, help ='The version of the genome')

    args = parser.parse_args()

    bed1 = os.path.abspath(args.bed1)
    
    
    try:
        os.path.isfile(bed1)
        f1=open(bed1,'r')
    except IOError:
        print('Error: There was no input whole genome variants file!')
        sys.exit()

 
    bed2 = os.path.abspath(args.bed2)
    
    try:
       os.path.isfile(bed2)
       f1=open(bed2,'r')
    except IOError:
        print('Error: There was no input the known brain eqtls file!')
        sys.exit()


    bPCA = os.path.abspath(args.bPCA)

    out = args.out
    hg=args.hg

   
    ##Converted the corrdnates to hg19

    def liftover():
        if (hg=="hg18") or (hg=="hg38"):
            liftOver(hg,bed1)
        else:
            pass
    liftover()


    print ( 'Searching for enhancers and promoters peaks...')
    get_annotate1(bed1, args.table_annovar_path)
    get_annotate2(bed1, args.table_annovar_path)
    get_annotate3(bed1, args.table_annovar_path)
    get_annotate4(bed1, args.table_annovar_path)

    os.system(args.table_annovar_path+'''/table_annovar.pl '''+bed1+' '+args.table_annovar_path+'''/humandb/ -buildver hg19 -out '''+out+''' -remove -protocol refGene,eigen,cadd,dann,gwava,fathmm,gnomad_genome -operation g,f,f,f,f,f,f -nastring .''')

    print "NOTICE: Start separate the exonic, intergenic,...,genes";
   
    ann_file= out+".hg19_multianno.txt"
    ann=open(ann_file, 'r')
    inter=open('intergenic_gene_regions.txt', 'w')
    up=open('updown_stream_gene_regions.txt', 'w')
    ex=open('exonic_gene_regions.txt', 'w')

    def sep():
        for ann1 in ann.readlines():
            ann2=ann1.strip()
            ann3=ann2.split('\t')
            if ann3[5]=="intergenic":
                ann33=ann3[0] +"\t" + ann3[1] +"\t" + ann3[2] +"\t" + ann3[3] +"\t" + ann3[4] +"\t" + ann3[5] +"\t" + ann3[6] +"\t" + ann3[7] +"\t" + ann3[8] +"\t" + ann3[9] +"\t" + ann3[10] +"\t" + ann3[11] +"\t" + ann3[12] +"\t" + ann3[13] +"\t" + ann3[14] +"\t" + ann3[15] +"\t" + ann3[16] +"\t" + ann3[17] +"\t" + ann3[18] +"\t" + ann3[19] +"\t" + ann3[20] +"\t" + ann3[21] +"\t" + ann3[22] +"\t" + ann3[23] +"\t" + ann3[24] +"\t" + ann3[25] +"\t" + ann3[26]
                inter.write('%s\n' %ann33)
            elif ann3[5]=="upstream;downstream":
                ann33=ann3[0] +"\t" + ann3[1] +"\t" + ann3[2] +"\t" + ann3[3] +"\t" + ann3[4] +"\t" + ann3[5] +"\t" + ann3[6] +"\t" + ann3[7] +"\t" + ann3[8] +"\t" + ann3[9] +"\t" + ann3[10] +"\t" + ann3[11] +"\t" + ann3[12] +"\t" + ann3[13] +"\t" + ann3[14] +"\t" + ann3[15] +"\t" + ann3[16] +"\t" + ann3[17] +"\t" + ann3[18] +"\t" + ann3[19] +"\t" + ann3[20] +"\t" + ann3[21] +"\t" + ann3[22] +"\t" + ann3[23] +"\t" + ann3[24] +"\t" + ann3[25] +"\t" + ann3[26]
                up.write('%s\n' %ann33)
            else:
                ann33=ann3[0] +"\t" + ann3[1] +"\t" + ann3[2] +"\t" + ann3[3] +"\t" + ann3[4] +"\t" + ann3[5] +"\t" + ann3[6] +"\t" + ann3[7] +"\t" + ann3[10] +"\t" + ann3[12] +"\t" + ann3[13] +"\t" + ann3[14] +"\t" + ann3[17] +"\t" + ann3[19] 
                ex.write('%s\n' %ann33)
        inter.close()
        up.close()
        ex.close()
    sep()
        


    

    f1 =open('intergenic_gene_regions.txt', 'r')

    gene_name=6 # The gene names column in annovar output file
    dist=7      # The distance column in annovar output file
    def intergenic(gene_name,dist):
        with open(out + "_intergenic_gene.txt", 'w') as output:
                for l1 in f1.readlines():
                    l2=l1.strip()
                    l3=l2.split()[gene_name]
                    l4=l3.split(',')
                    l5=l4[0]
                    d3=l2.split()[dist]
                    d4=d3.split(';')
                    d5=d4[0]
                    d6=d5.split('=')
                    d7=d6[1]
            
                    gene1=l2.split()[0]  + "\t" +  l2.split()[1]  + "\t" + l2.split()[2]  + "\t" + l2.split()[3]  + "\t" +  l2.split()[4]  + "\t" + l2.split()[5]  +  "\t" + l5 + "\t" +  d7 \
            + "\t" + l2.split()[10] + "\t" + l2.split()[12] + "\t" + l2.split()[13] + "\t" + l2.split()[14] + "\t" + l2.split()[17] + "\t" + l2.split()[19] 
            
                    g23=l1.split()[gene_name]
                    g24=g23.split(',')
                    g25=g24[1]
                    d23=l2.split()[dist]
                    d24=d23.split(';')
                    d25=d24[1]
                    d26=d25.split('=')
                    d27=d26[1]
                    gene2= l2.split()[0]  + "\t" +  l2.split()[1]  + "\t" + l2.split()[2]  + "\t" +  l2.split()[3]  + "\t" +  l2.split()[4]  + "\t" + l2.split()[5]  +  "\t" + g25 + "\t" +  d27 \
            + "\t" + l2.split()[10] + "\t" + l2.split()[12] + "\t" + l2.split()[13] + "\t" + l2.split()[14] + "\t" + l2.split()[17] + "\t" + l2.split()[19]       

                    output.write('%s\n' %gene1 )
                    output.write('%s\n' %gene2 )
                output.close()
    intergenic(gene_name,dist=7)

    

    up =open('updown_stream_gene_regions.txt', 'r')

    gene_name=6 # the gene names column in annovar output file
    up_down=5      # The distance column in annovar output file
   

    def up_down_stream(gene_name,up_down):
        with open(out + "_up_down_stream_gene.txt", 'w') as output:
            for l1 in up.readlines():
                l2=l1.strip()
                l3=l2.split()[gene_name]
                l4=l3.split(';')
                l5=l4[0]
                d3=l2.split()[up_down]
                d4=d3.split(';')
                d5=d4[0]
          
            

                gene1=l2.split()[0]  + "\t" +  l2.split()[1]  + "\t" + l2.split()[2]  + "\t" + l2.split()[3]  + "\t" +  l2.split()[4]   +  "\t" + d5 + "\t" +  l5 + "\t" + l2.split()[7] \
            + "\t" + l2.split()[10] + "\t" + l2.split()[12] + "\t" + l2.split()[13] + "\t" + l2.split()[14] + "\t" + l2.split()[17] + "\t" + l2.split()[19]

                g23=l1.split()[gene_name]
                g24=g23.split(';')
                g25=g24[1]
                d23=l2.split()[up_down]
                d24=d23.split(';')
                d25=d24[1]
                gene2= l2.split()[0]  + "\t" +  l2.split()[1]  + "\t" + l2.split()[2]  + "\t" +  l2.split()[3]  + "\t" +  l2.split()[4]   +  "\t" + d25 + "\t" +  g25 + "\t" + l2.split()[7]\
            + "\t" + l2.split()[10] + "\t" + l2.split()[12] + "\t" + l2.split()[13] + "\t" + l2.split()[14] + "\t" + l2.split()[17] + "\t" + l2.split()[19]

                output.write('%s\n' %gene1 )
                output.write('%s\n' %gene2 )
            output.close()
    
    up_down_stream(gene_name,up_down)
    
    os.system('''cat exonic_gene_regions.txt ''' +out + "_intergenic_gene.txt" +''' ''' +out + "_up_down_stream_gene.txt" +''' > Final_annovar_output_1.txt''')
    os.system('''sed '1d;2d' Final_annovar_output_1.txt > Final_annovar_output_id.txt''')
    os.system('''rm ''' +out + "_intergenic_gene.txt" +''' ''' +out + "_up_down_stream_gene.txt" +''' Final_annovar_output_1.txt ''')
     
    
    

    a=open('Final_annovar_output_id.txt','r')
    with open("variants_ids_gtex.txt", 'w') as outp:
        for j in a.readlines():
                j1=j.split('\t')
                j2=j1[0]
                if (j2=="chrX") or (j2=="chrY"):
                    j3=j2
                    j4=j3+"_"+j1[1]+"_"+j1[3]+"_"+ j1[4] +"_"+ hg19
                    outp.write('%s\n' %j4)
                else:
                    j5=j2
                    hg19='b37'
                    j6=j5+"_"+j1[1]+"_"+j1[3]+"_"+ j1[4] +"_"+ hg19
                    outp.write('%s\n' %j6)
        outp.close()
   
    os.system('''paste variants_ids_gtex.txt Final_annovar_output_id.txt > Final_annovar_output.txt''')
    
    


    print "NOTICE: Processing the known brain eQTLs score for each mutated gene\n";
    
    f2=open(bed2)
    dict_lines = {}
    
    gene_name=7
    phen_gene_name=2
    
    
    for p1 in f2.readlines():
        dict_lines[p1.split()[0]] = p1
    
    f3=open("Final_annovar_output.txt", 'r')
    
    with open(out +"_eqtls_score.txt", 'w') as output2:
        for val in f3.readlines():
            val1=val.strip()
            val2=val1.split()[0] 
            #val2=val1.split()[7]
            p1 = dict_lines.get(val2, None)
            if p1==None:
                p2= str(0.5)
                output2.write('%s\n' % p2)
            else:
                p3= p1.split()[2]
                output2.write('%s\n' % p3)
        output2.close()
 
   
#################################################
    
    print "NOTICE: Processing the enhancer/promtors: H3K27Ac";
    
    for file in glob.glob("H3K27Ac/*_bed.txt"):
        f = open(( file.rsplit( ".", 1 )[ 0 ] ) + "_bed", "w")
        file1=open(file,'r')
        f5=open("Final_annovar_output.txt", 'r')
        dict_lines = {}
        for g1 in file1.readlines():

            dict_lines[g1.split()[0]] = g1
        for g2 in f5.readlines():
            g3=g2.strip()
            g4=g3.split('\t')[0] 
            g5 = dict_lines.get(g4, None)
            if g5==None:
                g6=str(0)
                f.write('%s\n' % g6)
            else:
                g7= str(1)
                f.write('%s\n' % g7)
        f.close()


    os.system('paste H3K27Ac/*_bed_bed> H3K27Ac_score_1.txt')

    print "NOTICE: Processing the enhancer/promtors score H3K27me3";

  
    for file in glob.glob("H3K27me3/*_bed.txt"):
        f = open(( file.rsplit( ".", 1 )[ 0 ] ) + "_bed", "w")
        file1=open(file,'r')
        f5=open("Final_annovar_output.txt", 'r')
        dict_lines = {}
        for g1 in file1.readlines():
            
            dict_lines[g1.split()[0]] = g1
            
        for g2 in f5.readlines():
            g3=g2.strip()
            g4=g3.split('\t')[0]
            g5 = dict_lines.get(g4, None)
            if g5==None:
                g6=str(0)
                f.write('%s\n' % g6)
            else:
                g7= str(1)
                f.write('%s\n' % g7)
        f.close()

    
    os.system('paste H3K27me3/*_bed_bed > H3K27me3_score_1.txt')

    print "NOTICE: Processing the enhancer/promtors score for H3K4me1";

    for file in glob.glob("H3K4me1/*_bed.txt"):
        f = open(( file.rsplit( ".", 1 )[ 0] ) + "_bed", "w")
        file1=open(file,'r')
        f5=open("Final_annovar_output.txt", 'r')
        dict_lines = {}
        for g1 in file1.readlines():

            dict_lines[g1.split()[0]] = g1

        for g2 in f5.readlines():
            g3=g2.strip()
            g4=g3.split('\t')[0]
            g5 = dict_lines.get(g4, None)
            if g5==None:
                g6=str(0)
                f.write('%s\n' % g6)
            else:
                g7= str(1)
                f.write('%s\n' % g7)
        f.close()
  
    os.system('paste H3K4me1/*_bed_bed > H3K4me1_score_1.txt')

    ######################
    print "NOTICE: Processing the enhancer/promtors score for H3K4me3";

    for file in glob.glob("H3K4me3/*_bed.txt"):
        f = open(( file.rsplit( ".", 1 )[ 0] ) + "_bed", "w")
        file1=open(file,'r')
        f5=open("Final_annovar_output.txt", 'r')
        dict_lines = {}
        for g1 in file1.readlines():

            dict_lines[g1.split()[0]] = g1

        for g2 in f5.readlines():
            g3=g2.strip()
            g4=g3.split('\t')[0]
            g5 = dict_lines.get(g4, None)
            if g5==None:
                g6=str(0)
                f.write('%s\n' % g6)
            else:
                g7= str(1)
                f.write('%s\n' % g7)
        f.close()

    os.system('paste H3K4me3/*_bed_bed > H3K4me3_score_1.txt')

   ####################
 
  
    H3K27Ac_score_mean=open('H3K27Ac_score_1.txt', 'r')
    H3K27Ac_score_mean_out=open('H3K27Ac_score.txt', 'w')

    for s1 in H3K27Ac_score_mean:
        s2=s1.split()
        s3=np.array(s2).astype(np.float)
        s4=np.mean(s3)
        if s4 ==0:
            s5=s4
            H3K27Ac_score_mean_out.write('%s\n' % s5)
        else:
            s4=1
            H3K27Ac_score_mean_out.write('%s\n' % s4)
    H3K27Ac_score_mean_out.close()

    H3K27me3_score_mean=open('H3K27me3_score_1.txt', 'r')

    H3K27me3_score_mean_out=open('H3K27me3_score.txt', 'w')

    for s1 in H3K27me3_score_mean:
        s2=s1.split()
        s3=np.array(s2).astype(np.float)
        s4=np.mean(s3)
        if s4 ==0:
            s5=s4
            H3K27me3_score_mean_out.write('%s\n' % s5)
        else:
            s4=1
            H3K27me3_score_mean_out.write('%s\n' % s4)
    H3K27me3_score_mean_out.close()

    H3K4me1_score_mean=open('H3K4me1_score_1.txt', 'r')

    H3K4me1_score_mean_out=open('H3K4me1_score.txt', 'w')

    for s1 in H3K4me1_score_mean:
        s2=s1.split()
        s3=np.array(s2).astype(np.float)
        s4=np.mean(s3)
        if s4 ==0:
            s5=s4
            H3K4me1_score_mean_out.write('%s\n' % s5)
        else:
            s4=1
            H3K4me1_score_mean_out.write('%s\n' % s4)
    H3K4me1_score_mean_out.close()
 
    H3K4me3_score_mean=open('H3K4me3_score_1.txt', 'r')

    H3K4me3_score_mean_out=open('H3K4me3_score.txt', 'w')

    for s1 in H3K4me3_score_mean:
        s2=s1.split()
        s3=np.array(s2).astype(np.float)
        s4=np.mean(s3)
        if s4 ==0:
            s5=s4
            H3K4me3_score_mean_out.write('%s\n' % s5)
        else:
            s4=1
            H3K4me3_score_mean_out.write('%s\n' % s4)
    H3K4me3_score_mean_out.close()




    print "NOTICE: Processing the final variant score";

    filenames = ['''Final_annovar_output.txt''' , '''''' + out + "_eqtls_score.txt" +'''''','''H3K27Ac_score.txt''', '''H3K27me3_score.txt''', '''H3K27me3_score.txt''','''H3K4me3_score.txt''']

    with open('Final_variants_score.txt', 'w') as outfile:
        for fname in filenames:
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)


   # os.system('''paste Final_annovar_output.txt ''' + out + "_eqtls_score.txt" +''' H3K27Ac_score.txt H3K27me3_score.txt H3K4me1_score.txt H3K4me3_score.txt > Final_variants_score_1.txt''') 
    #os.system('''awk '{if ($13!=182) print $0}' Final_variants_score_1.txt > Final_variants_score.txt''')

    print "NOTICE: Processing the imputation";

    
    #subprocess.call ("/usr/bin/Rscript --vanilla imput_bpca_gene.R", shell=True)
    subprocess.call ("""Rscript --vanilla """+ bPCA + """  """, shell=True)

    print "NOTICE: Eigen missing values";
    os.system('''awk '{if ($10==".") print $0}' Final_variants_score.txt | wc -l''')

    print "NOTICE: CADD missing values";
    os.system('''awk '{if ($11==".") print $0}' Final_variants_score.txt | wc -l''')
    
    print "NOTICE: DANN missing values";
    os.system('''awk '{if ($12==".") print $0}' Final_variants_score.txt | wc -l''')    
    
    print "NOTICE: GWAVA missing values";
    os.system('''awk '{if ($13==".") print $0}' Final_variants_score.txt | wc -l''')

    print "NOTICE: FTAHMM missing values";
    os.system('''awk '{if ($14==".") print $0}' Final_variants_score.txt | wc -l''')

    print "NOTICE: Done; please check your results";

if __name__ == '__main__':
    main()






