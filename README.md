# ChromTransfer
Code for modeling, feature attribution analysis, and model interpretation as well as the resulting models used in Salvatore et al, 2022 (BioRxiv).

[![DOI](https://zenodo.org/badge/454701073.svg)](https://zenodo.org/badge/latestdoi/454701073)

Models were trained using tensorflow-gpu (version 1.14) and keras (version 3.2.1) on a Linux SMP Debian 4.19.208-1 x86_64 machine using NVIDIA Quadro RTX 6000 cards with 24 GB of VRAM.


## Finetuning

Finetuning may be performed on chromatin accessability data of new cell lines via the `finetune.py` script. 

```
usage: finetune.py [-h] --positive-sequence-name POSITIVE_SEQUENCE_NAME --pretrained-model <pretrained-model.h5> [--epochs EPOCHS] [--shuffle SHUFFLE] [--cache]
                   --output <model.h5>
                   <sequences.fasta>

positional arguments:
  <sequences.fasta>     Training sequences in FASTA format.

optional arguments:
  -h, --help            show this help message and exit
  --positive-sequence-name POSITIVE_SEQUENCE_NAME
                        Name of positive-class samples. E.g. 'A549' for headers of form '>A549::chr2:1527286-1527886'. Other names are treated as negatve instances.
  --pretrained-model <pretrained-model.h5>
                        Pretrained model.
  --epochs EPOCHS       Number of finetuning epochs.
  --shuffle SHUFFLE     Dataset shuffle buffer size.
  --cache               Cache dataset in memory.
  --output <model.h5>   Finetuned model destination.
```


## Using the models for predictions

We provide a stand-alone script that allows you to use the 7 models (pre-trained, A549, HCT116, HepG2, GM12878, K562, MCF7) for predictions of chromatin accessibility on a set of genomic regions (bed). 

First, you need to install the genomic_tool_tl.yml environment:
    
    conda env create -f genomic_tool_tl.yml 

The prediction script is run in the following way:
  
    python src/run_prediction.py -b data/tmp.bed -r genome_file/hg38seq.fa -g 0 -f tmp.fasta -o output_dir
  
where:
  
    -b = bed file (necessary)
    
    -r = genome file in fasta format (necessary)
    
    -g = make use of GPU(s) (1) or not (0)
    
    -f = file to store temporary fasta (necessary)
    
    -o = output folder to store csv result file (necessary)
    
