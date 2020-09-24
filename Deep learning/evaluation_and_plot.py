#%%
from pyteomics import mzml
import pandas as pd
import tensorflow as tf
import numpy as np
import time

# %%
max_len = 12 # maximum peptide length
k = 5 # k-mer length

pad = '_'
#alphabet = ['V', 'H', 'L', 'A', 'N', 'P', 'Y', 'W', 'T', 'E', 'Q', 'G', 'M', 'S', 'R', 'D', 'C', 'I', 'K', 'F',pad]
alphabet = [pad,'V', 'H', 'L', 'A', 'N', 'P', 'Y', 'W', 'T', 'E', 'Q', 'G', 'M', 'S', 'R', 'D', 'C', 'I', 'K', 'F']

### MODIFICATIONS #############################################
mod_aa_alphabet = ['O']
delta_masses = [15.99491463]
assert(len(set(alphabet).intersection(set(mod_aa_alphabet))) == 0)
massdelta_to_char = dict(zip(delta_masses,mod_aa_alphabet))
aa_alphabet = alphabet + mod_aa_alphabet
###############################################################
aa_index = range(len(aa_alphabet))
aa_to_index = dict(zip(aa_alphabet,aa_index))
index_to_aa = dict(zip(aa_index,aa_alphabet))

def to_index(seq,aa_to_index=aa_to_index):
    return np.array(list(map(lambda c: aa_to_index[c], list(seq))))

def one_hot(a, num_classes=len(aa_alphabet)):
    a = np.array(a)  
    return np.squeeze(np.eye(num_classes)[a.reshape(-1)])

def rev_one_hot(one_hot_vec):
    return np.argmax(one_hot_vec, axis=-1)

def to_aa(indices,index_to_aa=index_to_aa):
    return np.array(list(map(lambda c: index_to_aa[c], list(indices))))

def padding(seq,max_len=12,pad=pad):
  n = len(seq)
  if n < max_len:
    return seq+pad*(max_len-n)
  else:
    return seq[:max_len]

def kmers(seq, k):
  n = len(seq)
  #result = [] 
  #for i in range(n-k+1):
  #  result.append(seq[i:i+k])
  result = list(map(lambda i: seq[i:i+k],range(n-k+1)))
  return result

def integrate_modifications(seq,mod_list,massdelta_to_char=massdelta_to_char):
  '''
  Inputs: 
  seq -  pepitde sequence e.g. 'TKMPYMGMA'
  mod_lost -  list of tuples: (pos,mod_type)
    [(pos1, mass1), (pos2, mass2)] 
    [(3, 15.99491463), (8, 15.99491463)]

  Output:
  seq : peptide sequence with exchanged characters e.g. 'M'+delta_mass -> 'O'
  e.g. 'TKMPYMGMA' becomes 'TKOPYMGOA' 
  '''       
  for loc, mass in mod_list:
      seq = list(seq)
      seq[loc-1] = massdelta_to_char[mass]
      seq = ''.join(seq)
  return seq

def do_all_the_shit(seq,max_len=max_len,k=k):
  '''
  Do: padding+kmer+indices+onehot
  Input:
  sequence
  Output:
  array with shape: (n-k+1,k,alphabet_size)
  '''
  padded = padding(seq,max_len=max_len) # add padding
  if k is None:
    k = len(padded)
  kmer = kmers(padded,k) # create list of kmers
  kmer_indices= list(map(lambda char: to_index(char),kmer)) # transform list of kmers to indices
  kmer_one_hot = list(map(lambda index: one_hot(index),kmer_indices)) # one-hot encode
  return np.array(kmer_one_hot)

if __name__ == "__main__":
        
    seq = 'CGICGCGDTGEHHHEHEHEH' 
    sequences = ['CGICGCGDTGEHHHEHEHEH' ,'TKMPYMGMA']

    for seq in sequences:#map(lambda x: x,sequences): 
      max_len = 30
      k = 5            
      padded = padding(seq,max_len=max_len)
      kmer = kmers(padded,k)
      kmer_indices= list(map(lambda char: to_index(char),kmer))
      kmer_one_hot = list(map(lambda index: one_hot(index),kmer_indices))
      print(np.array(kmer_one_hot).shape)
      print(kmer_indices)

    ### Try encoding
    indices = to_index(seq,aa_to_index)
    one_hot_encoded = one_hot(indices,len(aa_alphabet))
    
    reversed_indices = rev_one_hot(one_hot_encoded)
    reversed_seq = to_aa(reversed_indices,index_to_aa)
    print(reversed_indices)
    print(indices)
    print(seq)
    print(''.join(reversed_seq))

    ##### Modifications: 
    print(integrate_modifications('TKMPYMGMA',[(3, 15.99491463), (8, 15.99491463)]))
# %%

AUTOTUNE = tf.data.experimental.AUTOTUNE


max_len=24

MZ_MAX=1900
SPECTRUM_RESOLUTION=2
k = 50

def set_k(new_k):
    global k
    k = new_k
    return k

def tf_preprocess_spectrum(mz,intensity):
    '''
    converts a peaks list (mz,intensity) into a dense spectrum.
    '''
    #global MZ_MAX, SPECTRUM_RESOLUTION
   
    n_spectrum = MZ_MAX * 10**SPECTRUM_RESOLUTION
    mz = mz*10**SPECTRUM_RESOLUTION
    
    # TODO: check this:
    indices = tf.math.floor(mz)
    indices = tf.cast(indices,tf.int64)


    uniq_indices, i = tf.unique(indices)
    # TODO: check what exactly to use here, sum, max, mean, ...
    uniq_values = tf.math.segment_max(intensity,i)

    # create as mask to truncate between min<mz<max
    # eliminate zeros:
    lower_bound = 100 * 10**SPECTRUM_RESOLUTION
    notzero_mask = tf.math.greater(uniq_indices,tf.zeros_like(uniq_indices)+lower_bound)    
    # truncate :
    trunc_mask = tf.math.less(uniq_indices,tf.zeros_like(uniq_indices)+n_spectrum)
    # put into joint mask:
    mask = tf.logical_and(notzero_mask,trunc_mask)
    # apply mask:
    uniq_indices = tf.boolean_mask(uniq_indices,mask)
    uniq_indices = uniq_indices - lower_bound
    uniq_values = tf.boolean_mask(uniq_values,mask)
    

    #### workaroud, cause tf.SparseTensor only works with tuple indices, so with stack zeros
    zeros = tf.zeros_like(uniq_indices)
    uniq_indices_tuples = tf.stack([uniq_indices, zeros],axis = 1)
    sparse_spectrum = tf.SparseTensor(indices = uniq_indices_tuples, values = uniq_values,dense_shape = [n_spectrum-lower_bound,1])
    dense_spectrum = tf.sparse.to_dense(sparse_spectrum)

    return dense_spectrum

# def normalize(intensities):
#     max_int = tf.reduce_max(intensities)
#     normalized = tf.log(intensities+1.)-tf.log(max_int+1)
#     return normalized

def tf_maxpool(dense):
    shape = dense.shape
    dense = tf.reshape(dense,[1,-1,1,1])
    #k = 100
    n_spectrum = int(shape[0])
    x, i = tf.nn.max_pool_with_argmax(dense,[1,k,1,1],[1,k,1,1],padding='SAME')
    i0 = tf.constant(np.arange(0,n_spectrum,k))
    i0 = tf.reshape(i0,[1,int(n_spectrum/k),1,1]) 
    i = i-i0
    return x,i

def tf_maxpool_with_argmax(dense,k):
    dense = tf.reshape(dense,[-1,k]) # K=columns
    x = tf.reduce_max(dense,axis=-1)
    i = tf.math.argmax(dense,axis=-1)
    return x,i

def ion_current_normalize(intensities):
    total_sum = tf.reduce_sum(intensities**2)
    normalized = intensities/total_sum
    return normalized

def parse(mz,intensity):
    '''
    converts a peaks list (mz,intensity) into a two-vector spectrum
    '''
    intensity = ion_current_normalize(intensity)
    
    spectrum_dense = tf_preprocess_spectrum(mz, intensity)
    
    x,i = tf_maxpool_with_argmax(spectrum_dense,k=k)
    x = tf.cast(x,tf.float32)
    i = tf.cast(i,tf.float32)
    i = i/tf.cast(k,tf.float32)
    spectrum_two_vec = tf.stack([x,i],axis=1)
    return  spectrum_two_vec
# %%
def fire_up_generator(file_path="./ptm21.hdf5",n=1,preshuffle=True):
    with pd.HDFStore(file_path) as hdf:
        keys = np.array(hdf.keys())
        df = pd.DataFrame()
        for key in keys:
            print(key)
            df = df.append(hdf.select(key=key))
    n = len(df)
    print('n datapoints:',len(df))
    if preshuffle:
      df = df.sample(frac=1)
      print('preshuffling done.')
    global index
    global j       
    index = 0
    j = 0

    def generator():    
        global index
        r_entry = df
        if index > (n-1):
          index=0
        seq = str(r_entry.iloc[index]['seq'])
        seq = np.array(do_all_the_shit(seq,max_len=max_len,k=None))        
        mz,i = np.array(r_entry.iloc[index]['mz']), np.array(r_entry.iloc[index]['intensities'])
        index+=1
        yield mz,i,seq
        
    
    def neg_generator():        
        # global index
        global j
        r_entry = df
        if index > (n-1):
          # index=0 #reset to zero
          j=0    #reset the index from the first element
        seq = str(r_entry.iloc[n-j-1]['seq'])   
        seq = np.array(do_all_the_shit(seq,max_len=max_len,k=None))        
        mz,i = np.array(r_entry.iloc[j]['mz']), np.array(r_entry.iloc[j]['intensities'])
        j+=1
        yield mz,i,seq
        
    return generator,neg_generator
# %%
def get_dataset(generators,batch_size=10,training=False):
    generator,neg_generator = generators 
    ds_pos = tf.data.Dataset.from_generator(generator,output_types=(tf.float32,tf.float32,tf.float32),output_shapes=(None,None,(None,max_len,22)))
    ds_neg = tf.data.Dataset.from_generator(neg_generator,output_types=(tf.float32,tf.float32,tf.float32),output_shapes=(None,None,(None,max_len,22)))
    ds_pos = ds_pos.map(lambda mz,intensities,seq: ((parse(mz,intensities),seq),1),num_parallel_calls=AUTOTUNE)
    ds_neg = ds_neg.map(lambda mz,intensities,seq: ((parse(mz,intensities),seq),0),num_parallel_calls=AUTOTUNE)
    ds = tf.data.experimental.sample_from_datasets([ds_pos,ds_neg])
    ds = ds.repeat()
    if training:
      ds = ds.shuffle(100000)
    ds = ds.batch(batch_size)
    return ds

# %%
# read data, select 10 pairs
x,indeces = next(iter(get_dataset(fire_up_generator("/home/ubuntu/data/jiahao/prova.hdf5"))))

#%%
# check which one is positive/negative
indeces

#%%
#  positive pair
pos_spec = x[0][0]

pos_seq = x[1][0]

#%%
# negative pair
neg_spec = x[0][1]

neg_seq = x[1][1]


#%%
# reload the saved models, change directory, reload seq_model and siam_model 
sequence_model = tf.keras.models.load_model('/home/ubuntu/data/jiahao/files/_model_relu/sequence_model',compile=False)
base_model = tf.keras.models.load_model('/home/ubuntu/data/jiahao/files/_model_relu/siamese_model',compile=False)


# define inputs, intermediate inputs and otputs
my_input = sequence_model.input
intermediate_out = sequence_model.output
base_out = base_model(intermediate_out)

# feature extraction, pass the matrices to the embedding model
model_seq = tf.keras.Model(inputs=my_input, outputs=base_out)


# %%
#load spectrum model and extract faetures

spectrum_model = tf.keras.models.load_model('/home/ubuntu/data/jiahao/files/_model_relu/spectrum_model',compile=False)
base_model = tf.keras.models.load_model('/home/ubuntu/data/jiahao/files/_model_relu/siamese_model',compile=False)

input_ = spectrum_model.input
intermediate_out_spec = spectrum_model.output
base_out_ = base_model(intermediate_out_spec)

model_spec = tf.keras.Model(inputs=input_, outputs=base_out_)

# %%
pos_seq = tf.reshape(pos_seq, [1,1,24,22])
pos_seq_features = model_seq.predict(pos_seq)

print(pos_seq_features)

# %%
pos_spec = tf.reshape(pos_spec,[1,3600,2])
pos_spec_features = model_spec.predict(pos_spec)

print(pos_spec_features)


# %%
import matplotlib.pyplot as plt
import matplotlib.cm as cm

plt.scatter(pos_seq_features,pos_spec_features)
plt.xlabel('sequence embedding')
plt.ylabel('spectrum embedding_3')
plt.grid(True)
#plt.savefig('pos_pair_3.png')
# %%
# euclidean distance
np.linalg.norm(pos_seq_features-pos_spec_features)


# %%
neg_seq = tf.reshape(neg_seq, [1,1,24,22])
neg_seq_features = model_seq.predict(neg_seq)

print(neg_seq_features)

# %%
neg_spec = tf.reshape(neg_spec,[1,3600,2])
neg_spec_features = model_spec.predict(neg_spec)

print(neg_spec_features)

# %%
import matplotlib.pyplot as plt
import matplotlib.cm as cm

plt.scatter(neg_seq_features,neg_spec_features)
plt.title('negative pair')
plt.xlabel('sequence embedding')
plt.ylabel('spectrum embedding')
plt.grid(True)
#plt.savefig('neg_pair.png')

# %%
np.linalg.norm(neg_seq_features-neg_spec_features)


"""
the below code is used to plot the technical replicates
"""


# %%
def get_dataset_1(generators,batch_size=1,training=False):
    generator,neg_generator = generators 
    ds_pos = tf.data.Dataset.from_generator(generator,output_types=(tf.float32,tf.float32,tf.float32),output_shapes=(None,None,(None,max_len,22)))
    ds_neg = tf.data.Dataset.from_generator(neg_generator,output_types=(tf.float32,tf.float32,tf.float32),output_shapes=(None,None,(None,max_len,22)))
    ds_pos = ds_pos.map(lambda mz,intensities,seq: ((parse(mz,intensities),seq),1),num_parallel_calls=AUTOTUNE)
    ds_neg = ds_neg.map(lambda mz,intensities,seq: ((parse(mz,intensities),seq),0),num_parallel_calls=AUTOTUNE)
    ds = tf.data.experimental.sample_from_datasets([ds_pos,ds_neg])
    ds = ds.repeat()
    if training:
      ds = ds.shuffle(100000)
    ds = ds.batch(batch_size)
    return ds

#%%
# first replicate
x_0,indeces_0 = next(iter(get_dataset_1(fire_up_generator("/home/ubuntu/data/jiahao/one_sequence_1.hdf5"))))

#%%
# second replicate
x_1,indeces_1 = next(iter(get_dataset_1(fire_up_generator("/home/ubuntu/data/jiahao/one_sequence_2.hdf5"))))

#%%
# third replicate
x_2,indeces_2 = next(iter(get_dataset_1(fire_up_generator("/home/ubuntu/data/jiahao/one_sequence_3.hdf5"))))
#%%
#  techinical replicates
pos_spec_0 = x[0][0]

pos_seq_0 = x[1][0]

pos_spec_1 = x_1[0][0]

pos_spec_2 = x_2[0][0]

# %%
pos_seq_0 = tf.reshape(pos_seq_0, [1,1,24,22])
pos_seq_features_0 = model_seq.predict(pos_seq_0)

print(pos_seq_features)

# %%
pos_spec_0 = tf.reshape(pos_spec_0,[1,3600,2])
pos_spec_features_0 = model_spec.predict(pos_spec_0)

print(pos_spec_features_0)

#%%
pos_spec_1 = tf.reshape(pos_spec_1,[1,3600,2])
pos_spec_features_1 = model_spec.predict(pos_spec_1)

print(pos_spec_features_1)

#%%
pos_spec_2 = tf.reshape(pos_spec_2,[1,3600,2])
pos_spec_features_2 = model_spec.predict(pos_spec_2)

print(pos_spec_features_2)

# %%
import matplotlib.pyplot as plt

plt.scatter(pos_seq_features_0, pos_spec_features_0, label=f'euclidean distance_1 = {0.351}', alpha=0.4, s=18, marker='D') 
plt.scatter(pos_seq_features_0, pos_spec_features_1, label=f'euclidean distance_2 = {0.328}', alpha=1, s=14, marker='^', c='red')
plt.scatter(pos_seq_features_0, pos_spec_features_2, label=f'euclidean distance_3 = {0.316}', alpha=0.7, s=10)
plt.xlabel('sequence embedding')
plt.ylabel('spectrum embeddings')
plt.grid(True)
plt.legend()
#plt.savefig('tech_replicates_plot.png')

