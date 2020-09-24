import numpy as np
import pandas as pd
import tensorflow as tf
import time

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
        #global index
        global j
        r_entry = df
        if j > (n-1):
          #index=0 #reset to zero
          j=0    #reset the index from the first element
        seq = str(r_entry.iloc[n-j-1]['seq'])   
        seq = np.array(do_all_the_shit(seq,max_len=max_len,k=None))        
        mz,i = np.array(r_entry.iloc[j]['mz']), np.array(r_entry.iloc[j]['intensities'])
        j+=1
        yield mz,i,seq
        
    return generator,neg_generator

def get_dataset(generators,batch_size=128,training=True):
    generator,neg_generator = generators 
    ds_pos = tf.data.Dataset.from_generator(generator,output_types=(tf.float32,tf.float32,tf.float32),output_shapes=(None,None,(None,max_len,22)))
    ds_neg = tf.data.Dataset.from_generator(neg_generator,output_types=(tf.float32,tf.float32,tf.float32),output_shapes=(None,None,(None,max_len,22)))
    ds_pos = ds_pos.map(lambda mz,intensities,seq: ((parse(mz,intensities),seq),1),num_parallel_calls=AUTOTUNE)
    ds_neg = ds_neg.map(lambda mz,intensities,seq: ((parse(mz,intensities),seq),0),num_parallel_calls=AUTOTUNE)
    ds = tf.data.experimental.sample_from_datasets([ds_pos,ds_neg])
    ds = ds.repeat()
    if training:
      ds = ds.shuffle(buffer_size=100000)
    ds = ds.batch(batch_size)
    ds = ds.repeat()
    return ds
