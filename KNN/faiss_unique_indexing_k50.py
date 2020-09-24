#%%
# Install the FAISS library with the python bindings
# Import the FAISS library and other dependencies
import faiss
import numpy as np
import pandas as pd

from datetime import datetime
#%%
# Predefine the dimension and size of the database
dimension = 32  

# Database : peptide embeddings
xb = np.load("/home/ubuntu/data/jiahao/new_sequence_embeddings.npy")

# Query : spectra embeddings
xq = np.load("/home/ubuntu/data/jiahao/new_spectra_embeddings.npy")

### Specify number of centroids/clusters
nlist = len(xq)                                      #Number of clusters  #175921 
#%%
print(xb.shape)
print(xq.shape)

#   Shape of unique peptide sequences : (6780082, 32) vector matrix
#   Shape of spectra : (175921, )
#   Note that all vector values are stored in the float 32 type


#%%
### Check available number of GPUs
ngpus = faiss.get_num_gpus()
print("number of GPUs:", ngpus)

#%%
### Build index
start_index_time = datetime.now()
quantiser = faiss.IndexFlatL2(dimension)                                       #Build CPU index
#index_flat = faiss.IndexFlatL2(dimension)                                    #Build CPU index
#index_flat = faiss.GpuIndexFlatL2(dimension)                                #Build GPU index
#index_gpu_flat = faiss.index_cpu_to_gpu(res, 0, index_cpu_flat)
#quantiser = faiss.index_cpu_to_all_gpus(index_flat)                          #Convert to GPU_index
index = faiss.IndexIVFFlat(quantiser, dimension, nlist, faiss.METRIC_L2)
end_index_time = datetime.now()
print('Build Index Duration: {}'.format(end_index_time - start_index_time))
#Build Index Duration: 0:00:00.289215
#%%
### Index is trained to create 'nlist' number of clusters
### View state of the index, whether is_trained training is completed
###  'ntotal' is the number of index data.
print(index.is_trained)                                             # False (some no training)
start_train_time = datetime.now()
index.train(xb)                                                     # train on the database vectors
end_train_time = datetime.now()
print(index.ntotal)                                                 # 0

print('Train Duration: {}'.format(end_train_time - start_train_time))
#Train Duration: 0:02:21.520537
#%%
### After the index is created, add and search operations can be carried out based on the index. 
### Add the vectors to the clusters
start_add_time = datetime.now()
index.add(xb)                                                      # add the vectors and update the index
end_add_time = datetime.now()
print(index.is_trained)                                            # True
print('Total number of vectors added to the index: ', index.ntotal)             # 6780082

print('Add_vector Duration: {}'.format(end_add_time - start_add_time))
#Add_vector Duration: 0:00:17.424560
# %%
### Perform a search on the index for the Query vectors.
# ‘k’ specifies the number of similar vectors to be returned from the visited clusters.

k = 120                                                 # Number of nearest neighbor(s) to search
start_search_time = datetime.now()
D, I = index.search(xq, k)                                # actual search against query
end_search_time = datetime.now()
print("FAISS is done")

print('Actual search Duration: {}'.format(end_search_time - start_search_time))
#Actual search Duration: 0:00:00.579271
#%%
### Visualize indexes and distances

print(I[:5])                                           # neighbors of the 5 first queries
print(I[-5:])                                         # neighbors of the 5 last queries
print('I: ', I.shape) 
#(175921, 50)

print(D[:5])
print(D[-5:])
print('D: ', D.shape) 
#(175921, 50)

# The search operation returns the row indexes in the vector store
##   of the k most similar vectors for each query vector
##   along with their respective distances.

#%%
### Read in unique peptides datafile
unique_dt = pd.read_hdf('/home/ubuntu/data/jiahao/unique_df.hdf5', key="df", mode="r")
unique_dt.head(10)
print(unique_dt.shape)
#(6780082, 2)


#%%
### Inverted index list i.e. Indexed unique peptides
indexed_dt = unique_dt.iloc[np.squeeze(I.flatten()),:]
print(indexed_dt.shape)
#(8796050, 2)
#%%
### Save file
indexed_dt.to_csv("final_indexed_uniquepeptides_k50.csv", sep="\t",index=False)
#%%
### Calculate hyperscores
hyp_score = 1/(D.flatten())
print(len(hyp_score))
#(8796050,)
# %%
### Read in Spectra datafile
spec_dat = pd.read_hdf('/home/ubuntu/data/jiahao/spectra_df.hdf5', key="df", mode="r")
print(spec_dat.shape)
#(175921, 3)

#%%
### Prepare mapping of spectra to peptide according to number of peptide neighbors
spec_dat = spec_dat.loc[spec_dat.index.repeat(k)]
spec_dat.head(5)
print(spec_dat.shape)
#(8796050, 3)
#%%
### Hardcode e-value as 0
evalue = [0.0000]*len(spec_dat)                               #8796050
print(len(evalue))
# %%
### Create final .tsv dataframe
final_knn_df = pd.DataFrame({
    '#SpecFile':list(spec_dat.iloc[:,0]),
    'Title':list(spec_dat.iloc[:,1]),
    'Peptide':list(indexed_dt.iloc[:,0]),
    'Protein':list(indexed_dt.iloc[:,1]),
    'Hyperscore':list(hyp_score),
    'Evalue': list(evalue)
    })
#%%
print(final_knn_df.shape)                                          #8796050 rows × 6 columns
final_knn_df.head(5)
print(type(final_knn_df.loc[0, 'Protein']))
print(type(final_knn_df.loc[0, 'Hyperscore']))

#%%
### Expand dataframe based on peptide-to-protein number
final_knn_df = final_knn_df.explode('Protein').reset_index(drop=True)

print(final_knn_df.shape)                                           #22915063 rows × 6 columns

#%%
### Transform the lists to strings
#final_knn_df['Protein'] = final_knn_df['Protein'].str.join(', ')
#final_knn_df['Hyperscore'] = pd.DataFrame([str(line).strip('[').strip(']') for line in final_knn_df['Hyperscore']])

# %%
#Save final dataframe
final_knn_df.to_csv("final_knn_result_k50.csv", index=False)


# %%
### Remove any non-typtic amino acid (Keep only peptides whose amino acid ends with 'K' or ' ')
knn_check_df = pd.read_csv("final_knn_result_k50.csv")
knn_check_df = knn_check_df[~knn_check_df['Peptide'].str.endswith('V', 'H')]
knn_check_df = knn_check_df[~knn_check_df['Peptide'].str.endswith('L', 'A')]
knn_check_df = knn_check_df[~knn_check_df['Peptide'].str.endswith('N', 'P')]
knn_check_df = knn_check_df[~knn_check_df['Peptide'].str.endswith('Y', 'W')]
knn_check_df = knn_check_df[~knn_check_df['Peptide'].str.endswith('T', 'E')]
knn_check_df = knn_check_df[~knn_check_df['Peptide'].str.endswith('Q', 'G')]
knn_check_df = knn_check_df[~knn_check_df['Peptide'].str.endswith('M')]
knn_check_df = knn_check_df[~knn_check_df['Peptide'].str.endswith('S', 'D')]
knn_check_df = knn_check_df[~knn_check_df['Peptide'].str.endswith('C', 'I')]
knn_check_df = knn_check_df[~knn_check_df['Peptide'].str.endswith('F')]
#'V', 'H', 'L', 'A', 'N', 'P', 'Y', 'W', 'T', 'E', 'Q', 'G', 'M', 'S', 'D', 'C', 'I', 'F'

#%%
knn_check_df['Peptide'] = knn_check_df['Peptide'].str.replace('I', 'J')
knn_check_df['Peptide'] = knn_check_df['Peptide'].str.replace('L', 'J')
print(knn_check_df.shape)
#(22736158, 6)
# %%
knn_check_df.to_csv("final_best_knn_k50.csv", index=False)
# %%
