# sslRNAD

![banner](https://github.com/hklz/sslRNAD/blob/main/img/sRNA.png)

[![Webserver](https://img.shields.io/badge/Webserver-available-blue)](http://www.kangzlab.cn/)
![License](https://img.shields.io/badge/License-MIT-orange)
[![DOI](https://zenodo.org/badge/501138098.svg)](https://zenodo.org/badge/latestdoi/501138098)


This repository is the original implementation of sslRNAD: The **S**ingle **S**tem **L**oop s**R**NA **D**esigner (WebSever: *ssl*-sRNA *de novo* designer).


**sslRNAD** is a program that allows users to design single stem loop sRNAs (*ssl*-sRNA) from scratch, to repress the expression of the gene of interest (target gene). This package contains an *ssl*-sRNA scoring function and related scripts to design automatically *ssl*-sRNA with predefined single stem loop structure and tunable regulatory activities. 

Numerous *ssl*-sRNAs can be designed by sslRNAD for a selected gene. sslRNAD also accepts batch target genes input (DNA sequence in fasta-format). In this situation, one *ssl*-sRNA with strong activity will be designed for each target gene. For each *ssl*-sRNA, sslRNAD design simultaneously two pairs of primers, which containing the DNA sequences of the *ssl*-sRNA, the interfaced promoter (customized) and regions homologous to the expression backbone (customized vector/plasmid). With the designed primers, users can do one-pot PCR to construct the expression vector of the *ssl*-sRNA. A simple demo is provided. 

## System Requirements

### OS and Software Dependencies

Ubuntu 18.04

Python v3.6.9 or above

[ViennaRNA](https://www.tbi.univie.ac.at/RNA/#download) v2.4.17 or above

[NUPACK](http://www.nupack.org/downloads) v4.0.0.23 or above

### Python modules
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

This version of the software has been tested on : Python v3.6.9 and PyCharm Community Edition 2020.2.3 is used for this to run the python scripts.

## Instructions

The original source code is in the [*methods*](https://github.com/hklz/sslRNAD/tree/main/methods) folder, and you can run sslRNAD on *cmd* directly using these scripts. No additional installation is needed if aforementioned dependencies are satisfied. The detail usages are demonstrated in the *Demo*  section. Make sure the ViennaRNA is set to permanent global environment variables and the *methods*  folder is in the same path of the demo scripts.

## Demos
Simple demos are provided ```/sslRNAD/Demo```.
The [*batch_input_test.fasta*](https://github.com/hklz/sslRNAD/blob/main/Demo/batch_input_test.fasta) contains a set of genes which can be imported to sslRNAD to create a *ssl*-sRNA library. The [*flank_sequence_test.fasta*](https://github.com/hklz/sslRNAD/blob/main/Demo/flank_sequence_test.fasta) contains the flank sequences of the *ssl*-sRNA-insertion-site on the selected vector and the selected interfaced promoter sequence, which are required for the design of primers. Experimental construction of the *ssl*-sRNA expression vectors can be performed via one-pot PCR with the designed primers and the selected vector (as template).

The results of *ssl*-sRNA library and related primers are in the [*Demo_result.xlsx*](https://github.com/hklz/sslRNAD/blob/main/Demo/Demo_results.xlsx). Please note different results will be returned from every running of the programs, because *ssl*-sRNA sequences are generated randomly.

Design of *ssl*-sRNA with desired activity for a single gene could be achieved by directly input the target gene sequence, repression level and the number of required *ssl*-sRNA in *cmd*.

```
cd ~/Demo/
#As for constructing a library of target genes, an excel file containing the custom *ssl*-sRNAs and related primers will be generated.
#Usage: python3 Demo_batch_design.py -s <flank_sequence.fasta> -i <batch_input_target_genes.fasta> -o <output_file_name>
$ python3 Demo_batch_design.py -s flank_sequence_test.fasta -i batch_input_test.fasta -i batch_input_test.fasta -o Demo_results

#As for creating *ssl*-sRNA with desired activity for a single gene, *ssl*-sRNA candidates will be directily print on the screen.
#Usage: python3 Programmable_strength_sRNA.py -i <Target_24_nt_sequence> -r <Repression_level> -t <Trials>
$ python3 Programmable_strength_sRNA.py -i ATGCAGTCATCGTAGCAGTCAGTC -r S -t 5
```
The demo results are as follows.


![image](https://github.com/hklz/sslRNAD/blob/main/img/Demo_output.png)


![image](https://github.com/hklz/sslRNAD/blob/main/img/Programmable_strength_sRNA_design.png)



## Online simplified version

sslRNAD is freely available at [http://www.kangzlab.cn/](http://www.kangzlab.cn/) as *ssl*-sRNA *de novo* designer.

