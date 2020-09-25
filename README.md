# Large Scale Metaproteomic Analysis powered by Machine Learning Metric Embedding Approach

This project consists of 3 submodules:
* deep learning
* knn similarity search
* taxIt/X!Tandem

***For detailed information, please see the README file in each folder***

### Deep learning
In this submodule, we develope a siamese convolutional network

### Knn similarity search
Here, we perform peptide identification using the Knnn similarity search

### TaxIt/X!Tandem
We use an existing search engine software and we compare the results with ours.
In particular, we make the following modifications to the configuration file (snakemake.config.json):
- threads number : 36
- memory : 16GB
- search engine : X!Tandem
- 'sample : Run1_U4_2000ng.mgf' , use this [link](https://www.ebi.ac.uk/pride/archive/projects/PXD006118) to download the file.
- 'database reference : small.fasta' , [link](https://owncloud.hpi.de/s/fa0aV3lp4Mu8Upq) for downloading.

HexBin Plot of KNN and Tandem Outputs.R is added for results visualization 
