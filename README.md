# transfer_learning_for_regulatory_genomics
Repo for transfer learning


## The method was tested on Linux SMP Debian 4.19.208-1 x86_64 with NVIDIA Quadro RTX 6000
    -  the models were trained using tensorflow-gpu==1.14 and keras

## Before start the prediction you need to install the tl_andersson.yml environment:
    
    - conda env create -f genomic_tool_tl.yml 

# To run the prediction you need to use the script run_prediction.py in this way:
  
  python src/run_prediction.py -b data/tmp.bed -r genome_file/hg38seq.fa -g 0 -f tmp.fasta -o output_dir
  
  where:
  
    -b = bed file (necessary)
    
    -r = genome file in fasta format (necessary)
    
    -g = if you have GPU(s) otherwise not include this flag
    
    -f = file to store temporary fasta (necessary)
    
    -o = output folder to store csv result file (necessary)
    
        
