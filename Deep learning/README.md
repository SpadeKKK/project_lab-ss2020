# Siamese convolutional neural network

## system requirements
Python > 3.6

## install dependencies
install wget, numpy, pandas, pyteomics, TensorFlow > 2.0

## Step 1 : download data

```
python download_files.py
```

## Step 2 : generate hdf5
split training data and validation data, merge mzid and mzml files to generate peptide-spectrum pairs in pandas dataframe, convert dataframes to hdf format

* create training data hdf5 file
```
python merge_training.py
```
* create validation data hdf5 file
```
python merge_test.py
```
