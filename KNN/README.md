# KNN similary Search

## Unique peptide list
Retrieve the list of unique peptides from the [small.fasta](https://owncloud.hpi.de/s/fa0aV3lp4Mu8Upq)
```
python filtered_unique_peptides.py
```

## Peptide embeddings
run following script to create the peptide embeddings as database for the knn similarity search, model stored in this [folder](https://github.com/jiahao95/proteomics/tree/master/Deep%20learning/_model_relu_32)
```
python sequence_embedder.py
```

## Knn similarity search
```
python faiss_unique_indexing_k50.py
```
