# sslRNAD

![banner](https://github.com/hklz/sslRNAD/blob/main/img/sRNA.png)

[![Webserver](https://img.shields.io/badge/Webserver-available-blue)](http://www.kangzlab.cn/)
![License](https://img.shields.io/badge/License-MIT-orange)


This repository is the original implementation of sslRNAD: The **S**ingle **S**tem **L**oop s**R**NA **D**esigner.


**sslRNAD** is a program that allows users to design single stem loop sRNAs (*ssl-sRNA*) from scratch, to repress the expression of the gene of interest (target gene). This package contains an *ssl-sRNA* scoring function and related scripts to design automatically *ssl-sRNA* with predefined single stem loop structure and tunable regulatory activities. 

Numerous *ssl-sRNAs* can be designed by sslRNAD for a selected gene. sslRNAD also accepts batch target genes input (DNA sequence in fasta-format). In this situation, one *ssl-sRNA* with strong activity will be designed for each target gene. For each *ssl-sRNA*, sslRNAD design simultaneously two pairs of primers, which containing the DNA sequences of the *ssl-sRNA*, the interfaced promoter (customized) and regions homologous to the expression backbone (customized vector/plasmid). With the designed primers, users can do one-pot PCR to construct the expression vector of the ssl-sRNA. A simple demo is provided. 

## System Requirements

### Software and OS Dependencies
```
Ubuntu 18.04
Python v3.6.9 or above
ViennaRNA v2.4.17 or above
NUPACK v4.0.0.23 or above
```

### Python module
```
biopython
openpyxl
getopt
primer3-py
NUPACK
subprocess
random
re
```

This version of software has been tested on : Python v3.6.9 and PyCharm Community Edition 2020.2.3 are used for this study, to run the python scripts.

## Usage and Demo
A simple demo is provided ```/sslRNAD/Demo```.
You can run the Demo scripts along with *methods folder* in the same path.

```
cd ~/Demo/
#As for constructing a library of target genes, an excel file containing the custom sRNAs and related primers will be generated.
#Usage: python3 Demo_batch_design.py -s <flank_sequence.fasta> -i <batch_input_target_genes.fasta> -o <output_file_name>
$ python3 Demo_batch_design.py -s flank_sequence_test.fasta -i batch_input_test.fasta -i batch_input_test.fasta -o Demo_results

#As for creating sRNA with desired activity for a single gene, artificial sRNA candidates will be directily print on the screen.
#Usage: python3 Programmable_strength_sRNA.py -i <Target_24_nt_sequence> -r <Repression_level> -t <Trials>
$ python3 Programmable_strength_sRNA.py -i ATGCAGTCATCGTAGCAGTCAGTC -r S -t 5
```
The demo results are as follows.


![image](https://github.com/hklz/sslRNAD/blob/main/img/Demo_output.png)


![image](https://github.com/hklz/sslRNAD/blob/main/img/Programmable_strength_sRNA_design.png)



## Online simplified version

sslRNAD is freely available at [http://www.kangzlab.cn/](http://www.kangzlab.cn/)

