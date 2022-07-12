# sslRNAD

![banner](https://github.com/hklz/sslRNAD/blob/main/img/sRNA.png)

[![license](https://img.shields.io/badge/License-Free-yellowgreen)](LICENSE)
[![Webserver](https://img.shields.io/badge/Webserver-available-blue)](http://www.kangzlab.cn/)

This repository is the original implementation of sslRNAD: The **S**ingle **S**tem **L**oop **RNA** **Designer** .
**sslRNAD** is a program that allowing users to rapidly construct custom sRNAs targeting at genes of interest. This package contains a ssl-sRNA scoring function based on multiple linear regression model and related scripts to automatically generate ssl-sRNA with desired structure and activity. Moreover, batch target genes in a fasta-format file are acceptable as inputs.


## Dependency

Ubuntu 18.04

Python v3.6.9 or above

ViennaRNA v2.4.17 or above

NUPACK v4.0.0.23 or above


## Auto-design pipeline


Target gene sequences are filtered and normalized to standard base-pairing boxse of ssl-sRNA in [*input_analysis.py*](https://github.com/hklz/sslRNAD/blob/main/input_anaylsis.py)

Every Base-pairing box will further connect to auto-designed custom scaffold to make functional ssl-sRNA in [*sRNA_generator.py*](https://github.com/hklz/sslRNAD/blob/main/sRNA_generator.py)


For easy construction of the ssl-sRNA library, desired ssl-sRNAs plasmids can be construted via two pairs of primer designed in [*primer_design_local.py*](https://github.com/hklz/sslRNAD/blob/main/sRNA_generator.py). (The flank sequences of the sRNA insert sites should be saved in the [*arms.fasta*](https://github.com/hklz/sslRNAD/blob/main/arms.fasta))


## Usage
A simple demo is provided ```/sslRNAD/exmaple```.

```
cd ~/sslRNAD/exmaple
python3 sRNA_generator arms.fasta test.fasta
```

## Citing this work



## Online simplified version

sslRNAD is freely available at [http://www.kangzlab.cn/](http://www.kangzlab.cn/)

