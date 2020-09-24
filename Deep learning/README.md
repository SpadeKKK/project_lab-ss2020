# Siamese convolutional neural network

## System requirements
Python > 3.6

## Install dependencies
install wget, numpy, pandas, pyteomics, TensorFlow > 2.0

## Step 1 : download data

```
python download_files.py
```

## Step 2 : generate hdf5
split data into training data and validation data, test data; merge mzid and mzml files to generate peptide-spectrum pairs in pandas dataframe, convert dataframes to hdf format

* create training data hdf5 file
```
python merge_train_set.py
```
* create validation data hdf5 file
```
python merge_validation_set.py
```

## Step 3 : tf-dataset
the script performs preprocessing (peptide one-hot encoding and peak list parsing), shuffle and batch in a tensorflow dataset
```
python tf_dataset.py
```

## Step 4 : training 
it contains the architecture, performs the training, and saves the model
```
python model.py
```

## Step 5 : prepare spectra for knn search
- follow this [link](https://www.ebi.ac.uk/pride/archive/projects/PXD006118) and download Run1_U4_2000ng.raw
- follow instructions in [link](https://github.com/compomics/ThermoRawFileParser) to install command line tool for converting .raw file into .mzml file
- run following command on terminal
```
mono ThermoRawFileParser.exe -i=/home/user_folder/Run1_U4_2000ng.raw -o=/home/user_folder/output_name/ -f=2
```
- 

