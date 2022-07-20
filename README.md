# sslRNAD

![banner](https://github.com/hklz/sslRNAD/blob/main/img/sRNA.png)

[![license](https://img.shields.io/badge/License-Free-yellowgreen)](LICENSE)
[![Webserver](https://img.shields.io/badge/Webserver-available-blue)](http://www.kangzlab.cn/)

This repository is the original implementation of sslRNAD: The **S**ingle **S**tem **L**oop s**R**NA **D**esigner.


**sslRNAD** is a program that allows users to design single stem loop sRNAs (*ssl-sRNA*) from scratch, to repress the expression of the gene of interest (target gene). This package contains an *ssl-sRNA* scoring function and related scripts to design automatically *ssl-sRNA* with predefined single stem loop structure and tunable regulatory activities. 

Numerous *ssl-sRNAs* can be designed by sslRNAD for a selected gene. sslRNAD also accepts batch target genes input (DNA sequence in fasta-format). In this situation, one *ssl-sRNA* with strong activity will be designed for each target gene. For each *ssl-sRNA*, sslRNAD design simultaneously two pairs of primers, which containing the DNA sequences of the *ssl-sRNA*, the interfaced promoter (customized) and regions homologous to the expression backbone (customized vector/plasmid). With the designed primers, users can do one-pot PCR to construct the expression vector of the ssl-sRNA. A simple demo is provided. 


## Dependency

Ubuntu 18.04

Python v3.6.9 or above

ViennaRNA v2.4.17 or above

NUPACK v4.0.0.23 or above


## Auto-design pipeline


Target gene sequences are filtered and normalized to standard base-pairing boxse of *ssl*-sRNA in [*input_analysis.py*](https://github.com/hklz/sslRNAD/blob/main/input_anaylsis.py)

Every Base-pairing box will further connect to auto-designed custom scaffold to make functional *ssl*-sRNA in [*sRNA_generator.py*](https://github.com/hklz/sslRNAD/blob/main/sRNA_generator.py)

For easy construction of the *ssl*-sRNA library, desired *ssl*-sRNAs plasmids can be construted via two pairs of primer designed in [*primer_design_local.py*](https://github.com/hklz/sslRNAD/blob/main/sRNA_generator.py). (The flank sequences of the sRNA insert sites should be saved in the [*arms.fasta*](https://github.com/hklz/sslRNAD/blob/main/arms.fasta))


## Usage
A simple demo is provided ```/sslRNAD/exmaple```.

```
cd ~/sslRNAD/exmaple
python3 sRNA_generator arms.fasta test.fasta
```

## Citing this work



## Online simplified version

sslRNAD is freely available at [http://www.kangzlab.cn/](http://www.kangzlab.cn/)

