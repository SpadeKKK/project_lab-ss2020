# Large Scale Metaproteomic Analysis powered by Machine Learning Metric Embedding Approach

This project consists of 3 submodules:
* deep learning
* knn similarity search
* taxIt/X!Tandem

***For detailed information, please see the README file in each folder***

### Deep learning
In this submodule, we develope a siamese convolutional network

### Knn similarity search
Here, we performe peptide identification using the Knnn similarity search

### TaxIt/X!Tandem
We use an existing search engine software and we compare the results with ours.
In particular, we make the following modifications to the configuration file (snakemake.config.json):
- threads number : 36
- memory : 16GB
- search engine : X!Tandem
- 'sample : run1 U4 2000ng.mgf' , use this [link](ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2017/05/PXD006118) to download.
- 'database reference : small.fasta' , [link](https://owncloud.hpi.de/s/fa0aV3lp4Mu8Upq) for downloading.
