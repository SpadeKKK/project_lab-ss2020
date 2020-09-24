#%%
import os
import urllib.request
from urllib.request import urlretrieve
import Bio
from Bio import SeqIO
import pyteomics
from pyteomics import fasta, parser
import pandas as pd
import tensorflow as tf
import numpy as np

#%%
# read unique peptide file
df = pd.read_hdf('/home/ubuntu/data/jiahao/unique_df.hdf5', key="df", mode="r")

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


#%%
uniq_list = df['peptide'].tolist()

# %%
encoded = []
for pep in uniq_list: 
    res = do_all_the_shit(pep, max_len=24, k=None)
    encoded.append(res)

#%%
# convert encoded into a numpy array
uniq_array = np.asarray(encoded)


# %%
# disable eager execution, make it faster
tf.compat.v1.disable_eager_execution()

#%%
# recall the model
sequence_model = tf.keras.models.load_model('/home/ubuntu/data/jiahao/files/_model_relu/sequence_model',compile=False)
base_model = tf.keras.models.load_model('/home/ubuntu/data/jiahao/files/_model_relu/siamese_model',compile=False)

#%%
my_input = sequence_model.input
intermediate_out = sequence_model.output
base_out = base_model(intermediate_out)
seq_model = tf.keras.Model(inputs=my_input, outputs=base_out)

#%%
# generate sequnce embeddings
predictions = seq_model.predict(uniq_array, batch_size=256)

# %%
# save embeddings
np.save('sequence_embeddings.npy', predictions)



