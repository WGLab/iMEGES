## iMEGES: integrated Mental-disorder Genome Score

iMEGES is software that prioritize the mental disease genes. 

## Introduction

# Brain variant score (Variant prioritization)  
The first step of the iMEGES take non-coding variants as input and prioritize the variants. We integrated five different scores from different predictors (such as EIGEN, CADD, DANN, GWAVA, FATHMM), GNOMAD frequency, known brain eQTLs from CommonMind and enhancer/promotors from PsychENCODE of the variants to prioritize the variants based on deep learning. 

# Disease specific brain gene score (Gene prioritization) 
The gene prioritization step of iMEGES mainly based on annotations for each variant from the first step of iMEGES variant prioritization and machine learning based brain score named as ncDeepBrain.
To characterize each gene, we annotated each of the variants with three genomic feature scores: 
Machine learning brain variant score, general score (RVIs, GTEx, haploinsufficient scores), Phenolyzerdisease specific score such as Phenolyzer, CNVsSS score and de novo mutation score. To link the non-coding variants with PhenolyzerPhenolyzer and RVIs score to a gene, we assigned each of these variants to it closet gene in term of genome distance, i.e., distance to gene is 100KB or less. Considering the fact that maybe some genes harbor more than one mutations, we consider all the mutations and prioritize each of the variants for the specific mental disorders genes. After integration of these pieces of information, the output of the second layer is the deep learning probability for each of the mutated gene, namely the iMEGES score, which measure the susceptibility potential for this gene.

This is the GitHub repository for the documentation of the iMEGES software, described in the paper listed below. If you like this repository, please click on the "Star" button on top of this page, to show appreciation to the repository maintainer. If you want to receive notifications on changes to this repository, please click the "Watch" button on top of this page.


## Dependency
Annovar http://annovar.openbioinformatics.org/en/latest/

Phenolyzer http://phenolyzer.wglab.org/



## Installation 

Please clone the repository into your computer:

    git clone https://github.com/WGLab/iMEGES

Then enter iMEGES directory:

    cd iMEGES
    
## Synopsis

    python ncDeepBrain.py --help


## License Agreement

By using the software, you acknowledge that you agree to the terms below:

For academic and non-profit use, you are free to fork, download, modify, distribute and use the software without restriction.

## Contact
Atlas Khan (ak4046@cumc.columbia.edu)

Kai Wang (kw2701@cumc.columbia.edu)

## Reference

Khan A, Wang K, **iMEGES: integrated Mental-disorder GEnome score for prioritizing the susceptibility genes in personal mental disorders**, In prepration.

## More information
Wang Genomics Lab Homepage (http://wglab.org/)



