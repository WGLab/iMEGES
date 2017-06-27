## iMEGES: integrated Mental-disorder Genome Score

iMEGES is software that prioritize the mental disease genes. 

## Introduction

The first layer of the iMEGES take non-coding variants as input and prioritize the variants in this layer. We integrated five different scores from different predictors and the frequency of the variants to prioritize the variants. The second layer of iMEGES prioritize the mutated genes output form the first layer. Based on the phenolyzer, GTEX and RVIs score, we will have the list of the mutated prioritize disease genes. After running iMEGES, we will have the list of the prioritize disease genes for mental specific disease, such as schizophrenia.


## Dependency
Annovar

Phenolyzer

## Installation 

Please clone the repository into your computer:

    git clone https://github.com/WGLab/iMEGES

Then enter iMEGES directory:

    cd iMEGES
    
## Synopsis

    python iMEGES_17.py --help


## License Agreement

By using the software, you acknowledge that you agree to the terms below:

For academic and non-profit use, you are free to fork, download, modify, distribute and use the software without restriction.


## Reference

Khan A, Wang K, **iMEGES: integrated Mental-disorder GEnome score for prioritizing the susceptibility genes in personal mental disorders**, In prepration.

## Contact
Atlas Khan (ak4046@cumc.columbia.edu)

Kai Wang (kw@cumc.columbia.edu)




