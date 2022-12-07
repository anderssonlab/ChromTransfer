# ChromTransfer
Code for modeling, feature attribution analysis, and model interpretation as well as the resulting models used in Salvatore et al, 2022 (BioRxiv).

Models were trained using tensorflow-gpu (version 1.14) and keras (version 3.2.1) on a Linux SMP Debian 4.19.208-1 x86_64 machine using NVIDIA Quadro RTX 6000 cards with 24 GB of VRAM.


## Using the pre-trained model for fine-tuning on a new dataset

We provide an example script (nn_architecture/transfer_learning_example.py) to illustrate how to fine-tune the pre-trained model on new data.

Steps to follow :

    -1: Load pre-trained model

    -2: Decide at which level you want to fine-tune the original architecture [see comment in the script in the function transfer_learning_model()]
    
    -3: Load the file with sequences and labels of interest
    
    -4: Run training/validation/test or k-fold cross validation


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
    
