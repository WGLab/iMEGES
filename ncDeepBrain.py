#!/usr/bin/env python

import sys, argparse, os
from itertools import izip
import glob
import re
import numpy as np
import subprocess

#import get_annotation_tracks_from_psychencode as get_annot

#from scripts import myIO

prog_name = 'iMEGES_v01_delta_TP_v02.py'



def main():
    parser = argparse.ArgumentParser(description='The intergenic genes', prog = prog_name)
    parser.add_argument('-f1', '--anno',  required = True, metavar = 'annovar file', type = str, help ='The annovar file')
    parser.add_argument('-o', '--out',  required = True, metavar = 'Define the name of output file', type = str, help ='The name of output file')
    parser.add_argument('-E', '--eqtls',  required = True, metavar = 'eqtls_list', type = str, help ='The list of eqtls')	
    #parser.add_argument('-R', '--RVis',  required = True, metavar = 'RVis-score', type = str, help ='GTex_score_all.txt')    
   # parser.add_argument('-G', '--GTex',  required = True, metavar = 'GTex_score_all.txt', type = str, help ='GTex_score_all.txt')
    #parser.add_argument('-f2', '--anno',  required = True, metavar = 'annovar file', type = str, help ='The annovar file')
    args = parser.parse_args()
    anno = os.path.abspath(args.anno)
    out = os.path.abspath(args.out)
    #anno = os.path.abspath(args.anno)

    eqtls = os.path.abspath(args.eqtls)
  
    #get_annot.get_annotate(anno)

    os.system('''table_annovar.pl '''+ anno +'''  humandb/ -buildver hg19 -out '''+out+''' -remove -protocol refGene,eigen,cadd,dann,gwava,fathmm,gnomad_genome -operation g,f,f,f,f,f,f -nastring .''')

    print "NOTICE: Start separate the exonic, intergenic,...,genes";
    os.system('''awk '{if ($6=="intergenic" ) print $0}' ''' +out + ".hg19_multianno.txt" +'''   > intergenic_gene_regions.txt ''')
    os.system('''awk '{if ($6=="upstream;downstream" ) print $0}' ''' +out + ".hg19_multianno.txt" +'''  > updown_stream_gene_regions.txt ''')
    os.system('''awk '{if ($6!="intergenic" && $6!= "upstream;downstream") print $1 "\t"  $2  "\t" $3 "\t" $4 "\t" $5  "\t" $6  "\t" $7  "\t" $8 "\t" $11 "\t" $12 "\t" $14 \
    "\t" $15  "\t" $18 "\t" $20}' ''' +out + ".hg19_multianno.txt" +'''  > exonic_gene_regions.txt ''')
   

    

    f1 =open('intergenic_gene_regions.txt', 'r')

    gene_name=6 # the gene names column in annovar output file
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
            + "\t" + l2.split()[11] + "\t" + l2.split()[12] + "\t" + l2.split()[14] + "\t" + l2.split()[15] + "\t" + l2.split()[18] + "\t" + l2.split()[20] 
            
            	    g23=l1.split()[gene_name]
            	    g24=g23.split(',')
            	    g25=g24[1]
            	    d23=l2.split()[dist]
            	    d24=d23.split(';')
            	    d25=d24[1]
            	    d26=d25.split('=')
            	    d27=d26[1]
            	    gene2= l2.split()[0]  + "\t" +  l2.split()[1]  + "\t" + l2.split()[2]  + "\t" +  l2.split()[3]  + "\t" +  l2.split()[4]  + "\t" + l2.split()[5]  +  "\t" + g25 + "\t" +  d27 \
            + "\t" + l2.split()[11] + "\t" + l2.split()[12] + "\t" + l2.split()[14] + "\t" + l2.split()[15] + "\t" + l2.split()[18] + "\t" + l2.split()[20]       

            	    output.write('%s\n' %gene1 )
            	    output.write('%s\n' %gene2 )
        	output.close()
    intergenic(gene_name,dist=7)

    

    up =open('updown_stream_gene_regions.txt', 'r')

    gene_name=6 # the gene names column in annovar output file
    up_down=5      # The distance column in annovar output file
   

    def up_down_stream(gene_name,dist):
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
            + "\t" + l2.split()[11] + "\t" + l2.split()[12] + "\t" + l2.split()[14] + "\t" + l2.split()[15] + "\t" + l2.split()[18] + "\t" + l2.split()[20]

            	g23=l1.split()[gene_name]
            	g24=g23.split(';')
            	g25=g24[1]
            	d23=l2.split()[up_down]
            	d24=d23.split(';')
            	d25=d24[1]
            	gene2= l2.split()[0]  + "\t" +  l2.split()[1]  + "\t" + l2.split()[2]  + "\t" +  l2.split()[3]  + "\t" +  l2.split()[4]   +  "\t" + d25 + "\t" +  g25 + "\t" + l2.split()[7]\
            + "\t" + l2.split()[11] + "\t" + l2.split()[12] + "\t" + l2.split()[14] + "\t" + l2.split()[15] + "\t" + l2.split()[18] + "\t" + l2.split()[20]

            	output.write('%s\n' %gene1 )
            	output.write('%s\n' %gene2 )
            output.close()
    
    up_down_stream(gene_name,dist)
    
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
    
    f2=open(eqtls)
    dict_lines = {}
    
    gene_name=7
    phen_gene_name=2
   


    # def phenolyzer_score():
    for p1 in f2.readlines():
        dict_lines[p1.split()[phen_gene_name]] = p1

    f3=open("Final_annovar_output.txt", 'r')

    with open(out +"_eqtls_score.txt", 'w') as output2:
        for val in f3.readlines():
            val1=val.strip()
	    val2=val1.split()[7]
            p1 = dict_lines.get(val2, None)
            if p1==None:
                p2= str(0.5)
                output2.write('%s\n' % p2)
            else:
                p3= p1.split()[3]
                output2.write('%s\n' % p3)
        output2.close()

    
    
    

#################################################


    #get_annot.get_annotate(anno)

    
    print "NOTICE: Processing the enhancer/promtors score for each mutated gene in CNON";
    
    #H3K27Ac="H3K27Ac"
    #os.chdir(H3K27Ac)

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

    


    #os.chdir("../")
    os.system('paste H3K27Ac/*_bed_bed> H3K27Ac_score_1.txt')


    print "NOTICE: Processing the enhancer/promtors score for each mutated gene in EpiMap";

    #H3K27me3="H3K27me3"
    #os.chdir(H3K27me3)
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

    #os.chdir("../")
    os.system('paste H3K27me3/*_bed_bed > H3K27me3_score_1.txt')

    print "NOTICE: Processing the enhancer/promtors score for each mutated gene in Yale-ASD";

    #H3K4me1="H3K4me1"
    #os.chdir(H3K4me1)

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
  
    #os.chdir("../")
    os.system('paste H3K4me1/*_bed_bed > H3K4me1_score_1.txt')



    ######################
    print "NOTICE: Processing the enhancer/promtors score for each mutated gene in Yale-ASD";

    #H3K4me3="H3K4me3"
    #os.chdir(H3K4me3)

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

    #os.chdir("../")
    os.system('paste H3K4me3/*_bed_bed > H3K4me3_score_1.txt')

   ####################
 
  
    GTEX_score_mean=open('H3K27Ac_score_1.txt', 'r')
    GTEX_score_mean_out=open('H3K27Ac_score.txt', 'w')

    for s1 in GTEX_score_mean:
        s2=s1.split()
        s3=np.array(s2).astype(np.float)
        s4=np.mean(s3)
        if s4 ==0:
            s5=s4
            GTEX_score_mean_out.write('%s\n' % s5)
        else:
            s4=1
            GTEX_score_mean_out.write('%s\n' % s4)
    GTEX_score_mean_out.close()

    GTEX_score_mean=open('H3K27me3_score_1.txt', 'r')

    GTEX_score_mean_out=open('H3K27me3_score.txt', 'w')

    for s1 in GTEX_score_mean:
        s2=s1.split()
        s3=np.array(s2).astype(np.float)
        s4=np.mean(s3)
        if s4 ==0:
            s5=s4
            GTEX_score_mean_out.write('%s\n' % s5)
        else:
            s4=1
            GTEX_score_mean_out.write('%s\n' % s4)
    GTEX_score_mean_out.close()

    GTEX_score_mean=open('H3K4me1_score_1.txt', 'r')

    GTEX_score_mean_out=open('H3K4me1_score.txt', 'w')

    for s1 in GTEX_score_mean:
        s2=s1.split()
        s3=np.array(s2).astype(np.float)
        s4=np.mean(s3)
        if s4 ==0:
            s5=s4
            GTEX_score_mean_out.write('%s\n' % s5)
        else:
            s4=1
            GTEX_score_mean_out.write('%s\n' % s4)
    GTEX_score_mean_out.close()


 
    GTEX_score_mean=open('H3K4me3_score_1.txt', 'r')

    GTEX_score_mean_out=open('H3K4me3_score.txt', 'w')

    for s1 in GTEX_score_mean:
        s2=s1.split()
        s3=np.array(s2).astype(np.float)
        s4=np.mean(s3)
        if s4 ==0:
            s5=s4
            GTEX_score_mean_out.write('%s\n' % s5)
        else:
            s4=1
            GTEX_score_mean_out.write('%s\n' % s4)
    GTEX_score_mean_out.close()




    print "NOTICE: Processing the final variant score";

    os.system('''paste Final_annovar_output.txt ''' + out + "_eqtls_score.txt" +''' H3K27Ac_score.txt H3K27me3_score.txt H3K4me1_score.txt H3K4me3_score.txt > Final_variants_score_1.txt''') 
    os.system('''awk '{if ($14!=182) print $0}' Final_variants_score_1.txt > Final_variants_score.txt''')

    print "NOTICE: Processing the imputation";

    subprocess.call ("/usr/bin/Rscript --vanilla imput_bpca_gene.R", shell=True)
    print "NOTICE: Done; please check your results";



if __name__ == '__main__':
    main()






