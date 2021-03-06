## iMEGES: integrated Mental-disorder Genome Score

iMEGES is a software tool that prioritize the whole genome variants and genes related to mental disease. 

## Introduction
iMEGES has main two steps, which are as follows:

# Step 1: Brain variant score (Variant prioritization)  
The first step of the iMEGES take non-coding variants as input and prioritize the brain non-coding variants. We integrated five different scores from different predictors (such as EIGEN, CADD, DANN, GWAVA, FATHMM), GNOMAD frequency, known brain eQTLs from CommonMind and enhancer/promoters from PsychENCODE and Roadmap Epigenomics projects to prioritize the non-coding variants based on deep learning. 

# Step 2: Brain gene score (Gene prioritization) 
The gene prioritization step 2 of iMEGES mainly based on annotations for each variant from the first step of iMEGES; whole genome variant prioritization, brain score named as ncDeepBrain. To characterize each gene, we annotated each of the variants with three genomic feature scores: 
Machine learning brain variant score ncDeepBrain, general score (RVIs, GTEx, haploinsufficient scores), disease specific score such as Phenolyzer, CNVs score and denovo mutation score. To link the non-coding variants to a gene, we assigned each of these variants to it closet gene in term of genome distance, i.e., distance to gene is 100 KB or less. Considering the fact that maybe some genes harbor more than one mutations, we consider all the mutations and prioritize each of the variants for the specific mental disorders genes. After integration of these pieces of information, the output of the second Step is the deep learning probability for each of the mutated gene and call it the iMEGES gene score, which measure the susceptibility potential for this gene to be involved in mental disorders.

This is the GitHub repository for the documentation of the iMEGES software, described in the paper listed below. If you like this repository, please click on the "Star" button on top of this page, to show appreciation to the repository maintainer. If you want to receive notifications on changes to this repository, please click the "Watch" button on top of this page.


## Dependency

Annovar http://annovar.openbioinformatics.org/en/latest/

Phenolyzer http://phenolyzer.wglab.org/

EIGEN http://www.columbia.edu/~ii2135/eigen.html

CADD http://cadd.gs.washington.edu/

DANN https://cbcl.ics.uci.edu/public_data/DANN/ 

GWAVA https://www.sanger.ac.uk/sanger/StatGen_Gwava

FATHMM http://fathmm.biocompute.org.uk/

GNOMAD http://gnomad.broadinstitute.org/

PsychEncode https://www.synapse.org//#!Synapse:syn4921369/wiki/235539

EpiMap https://www.synapse.org/#!Synapse:syn5691268

CNON https://www.synapse.org/#!Synapse:syn4598822

Yale-ASD https://www.synapse.org/#!Synapse:syn4566311

Brain eQTLs https://www.synapse.org/#!Synapse:syn4622659

RVIs http://genic-intolerance.org/

GTEx http://www.gtexportal.org/home/

## Links of useful database

SZGR: http://bioinfo.mc.vanderbilt.edu/SZGR/

SZGene: http://www.szgene.org/

BDgene: http://bdgene.psych.ac.cn/index.do

BrainCloud: http://braincloud.jhmi.edu/

BrainSpan: http://www.brainspan.org/

NPdenovo: http://www.wzgenomics.cn/NPdenovo/

SZBD: https://academic.oup.com/schizophreniabulletin/article/43/2/459/2503611/SZDB-A-Database-for-Schizophrenia-Genetic-Research

SFARI GENE: https://gene.sfari.org/

AutDB: http://autism.mindspec.org/autdb/Welcome.do

AGD: http://autism-genetic-db.net/site/home 

AutismKB: http://autismkb.cbi.pku.edu.cn/index.php

International Genomics of Alzheimer's Project (IGAP): http://web.pasteur-lille.fr/en/recherche/u744/igap/igap_download.php

Roadmap Epigenomics: Brain enhancers and promoters
ftp://ftp.ncbi.nlm.nih.gov/pub/geo/DATA/roadmapepigenomics/by_sample/


## Computations related links

Keras: https://keras.io/

Tensorflow: https://www.tensorflow.org/

Theano: http://deeplearning.net/software/theano/

Scikit: http://scikit-learn.org/stable/documentation.html

## Tensorflow

We used keras library with Tensorflow backend for deep learning in our study. 

## Installation 

Please clone the repository into your computer:

    git clone https://github.com/WGLab/iMEGES

Then enter iMEGES directory:

    cd iMEGES
    
## Synopsis

## OPTIONS

* -h, --help show this help message and exit

*  -f1 input.bed, --bed1 input.bed
                        The Annovar input/bed format file
* -E input.bed, --bed2 input.bed
                        The list of the known brain eqtls from CommonMind
*  -B bPCA.R, --bPCA bPCA.R
                        The bPCA R code to do the imputation
*  -o The name of output file, --out The name of output file
                        The name of the output file
*  -buildver out_hg_name, --hg out_hg_name
                        The version of the genome


           python ncDeepBrain.py --help


## License Agreement

By using the software, you acknowledge that you agree to the terms below:

For academic and non-profit use, you are free to fork, download, modify, distribute and use the software without restriction.

## Contact
Atlas Khan (atlas.akhan@gmail.com)

Kai Wang (kaichop@gmail.com)

## Reference
Khan A, Wang K, **A deep learning based scoring system for prioritizing susceptibility variants for mental disorders**, IEEE International Conference on Bioinformatics and Biomedicine (BIBM). Page: 1698-1705, 2017, DOI: 10.1109/BIBM.2017.8217916

Khan A, Liu Q, Wang K, **iMEGES: integrated Mental-disorder GEnome score for prioritizing the susceptibility genes in personal genomes**, ACCEPTED BMC BIOINFOMARTICS.

## More information
Wang Genomics Lab Homepage (http://wglab.org/)



