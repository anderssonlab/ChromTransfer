# ChromTransfer
Code for modeling, feature attribution analysis, and model interpretation as well as the resulting models used in Salvatore et al, 2022 (BioRxiv).

Models were trained using tensorflow-gpu==1.14 and keras (version 3.2.1) on a Linux SMP Debian 4.19.208-1 x86_64 machine using NVIDIA Quadro RTX 6000 cards with 24 GB of VRAM.
    
## Using the models for predictions

We provide stand-alone script for using the models for predictions. 

First, you need to install the genomic_tool_tl.yml environment:
    
    - conda env create -f genomic_tool_tl.yml 

The prediction script is run in the following way:
  
  python src/run_prediction.py -b data/tmp.bed -r genome_file/hg38seq.fa -g 0 -f tmp.fasta -o output_dir
  
  where:
  
    -b = bed file (necessary)
    
    -r = genome file in fasta format (necessary)
    
    -g = to make use of GPU(s) (1) or not (0)
    
    -f = file to store temporary fasta (necessary)
    
    -o = output folder to store csv result file (necessary)
