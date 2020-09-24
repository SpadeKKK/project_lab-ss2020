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
### Check model performance
Once the training process is completed, run the following script to compute euclidean distance and generate plots
```
python evaluation_and_plot.py
```


## Step 5 : prepare spectra for the knn search
- follow this [link](https://www.ebi.ac.uk/pride/archive/projects/PXD006118) and download Run1_U4_2000ng.raw
- follow instructions in [link](https://github.com/compomics/ThermoRawFileParser) to install command line tool for converting .raw file into .mzml file
- run following command on terminal
```
mono ThermoRawFileParser.exe -i=/home/user_folder/Run1_U4_2000ng.raw -o=/home/user_folder/output_name/ -f=2
```
- run following python script to convert the .mzml file in hdf format
```
python convert_raw_spectra.py
```

## Step 6 : spectra embeddings
the following python script generates the spectra embeddings for the knn search, the model is stored in this [folder](https://github.com/jiahao95/proteomics/tree/master/Deep%20learning/_model_relu_32)
```
python spectra_embedder.py
```

***Continue with creating peptide embeddings and knn search in [here](https://github.com/jiahao95/proteomics/tree/master/KNN)***
