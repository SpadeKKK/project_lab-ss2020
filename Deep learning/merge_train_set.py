#%%
import os
from pyteomics import mzid, mzml
import pandas as pd
import numpy as np
import glob

"""
Files are downloaded and manually randomly divided into different folders
the following code is repeated but has the same effect, it is applied to various folders to 
generate pandas data frames and to store all the data in a single hdf5 file 
"""
#%%
os.chdir('./files/train')

mzid_files=glob.glob('*.mzid')
indexed_mzid = mzid.chain.from_iterable(mzid_files, use_index=True)

def _parse_mzid_entry(entry):
    spectrum_id = str(entry['spectrumID'])
    seq = str(entry['SpectrumIdentificationItem'][0]['PeptideSequence'])
    try: 
        mods = entry['SpectrumIdentificationItem'][0]['Modification'] 
    except:
        mods = None
    rank = int(entry['SpectrumIdentificationItem'][0]['rank'])
    file_location = str(entry['name'])
    return file_location,spectrum_id,seq,mods,rank

all_mzid = []
for entry in(indexed_mzid):
    all_mzid.append(_parse_mzid_entry(entry))
  

file_location,spectrum_ids,seq,mods,rank = zip(*all_mzid)

mzid_df = pd.DataFrame({'file':file_location,'id':spectrum_ids,'seq':seq})

def _parse_mzml_entry(entry):
    ID = str(entry['id'])
    mz = np.array(entry['m/z array'])
    intensities = np.array(entry['intensity array'])
    return ID, mz, intensities

all_spectra = []

for file in np.unique(file_location):
    print(file)
    indexed = mzml.MzML(file)
    for i,entry in enumerate(indexed.map(_parse_mzml_entry)):
        tupl = (file,)+entry
        all_spectra.append(tupl)

mzml_location, ids, mz, intensities = zip(*all_spectra)

spectra_df = pd.DataFrame({'file':mzml_location,'id':ids,'mz':mz,'intensities':intensities})

#### MERGE: mzid + mzml

merged_df = pd.merge(mzid_df,spectra_df,how='left',on=['file','id'])

merged_df = merged_df[['id','seq','mz','intensities']]

#%%
hdf = pd.HDFStore('/home/ubuntu/data/jiahao/files/train.hdf5', mode="w")

hdf.put(value=merged_df, key="df")

#%%
os.chdir('./train_1')

mzid_files_1=glob.glob('*.mzid')
indexed_mzid_1 = mzid.chain.from_iterable(mzid_files_1, use_index=True)


def _parse_mzid_entry(entry):
    spectrum_id = str(entry['spectrumID'])
    seq = str(entry['SpectrumIdentificationItem'][0]['PeptideSequence'])
    try: 
        mods = entry['SpectrumIdentificationItem'][0]['Modification'] 
    except:
        mods = None
    rank = int(entry['SpectrumIdentificationItem'][0]['rank'])
    file_location = str(entry['name'])
    return file_location,spectrum_id,seq,mods,rank

all_mzid_1 = []
for entry in(indexed_mzid_1):
    all_mzid_1.append(_parse_mzid_entry(entry))
  

file_location,spectrum_ids,seq,mods,rank = zip(*all_mzid_1)

mzid_df_1 = pd.DataFrame({'file':file_location,'id':spectrum_ids,'seq':seq})

def _parse_mzml_entry(entry):
    ID = str(entry['id'])
    mz = np.array(entry['m/z array'])
    intensities = np.array(entry['intensity array'])
    return ID, mz, intensities

all_spectra_1 = []

for file in np.unique(file_location):
    print(file)
    indexed = mzml.MzML(file)
    for i,entry in enumerate(indexed.map(_parse_mzml_entry)):
        tupl = (file,)+entry
        all_spectra_1.append(tupl)

mzml_location, ids, mz, intensities = zip(*all_spectra_1)

spectra_df_1 = pd.DataFrame({'file':mzml_location,'id':ids,'mz':mz,'intensities':intensities})

#### MERGE: mzid + mzml

merged_df_1 = pd.merge(mzid_df_1,spectra_df_1,how='left',on=['file','id'])

merged_df_1 = merged_df_1[['id','seq','mz','intensities']]

#%%
hdf.put(value=merged_df_1, key="df1")

# %%
os.chdir('./train_2')

mzid_files_2=glob.glob('*.mzid')
indexed_mzid_2 = mzid.chain.from_iterable(mzid_files_2, use_index=True)


def _parse_mzid_entry(entry):
    spectrum_id = str(entry['spectrumID'])
    seq = str(entry['SpectrumIdentificationItem'][0]['PeptideSequence'])
    try: 
        mods = entry['SpectrumIdentificationItem'][0]['Modification'] 
    except:
        mods = None
    rank = int(entry['SpectrumIdentificationItem'][0]['rank'])
    file_location = str(entry['name'])
    return file_location,spectrum_id,seq,mods,rank

all_mzid_2 = []
for entry in(indexed_mzid_2):
    all_mzid_2.append(_parse_mzid_entry(entry))
  

file_location,spectrum_ids,seq,mods,rank = zip(*all_mzid_2)

mzid_df_2 = pd.DataFrame({'file':file_location,'id':spectrum_ids,'seq':seq})

def _parse_mzml_entry(entry):
    ID = str(entry['id'])
    mz = np.array(entry['m/z array'])
    intensities = np.array(entry['intensity array'])
    return ID, mz, intensities

all_spectra_2 = []

for file in np.unique(file_location):
    print(file)
    indexed = mzml.MzML(file)
    for i,entry in enumerate(indexed.map(_parse_mzml_entry)):
        tupl = (file,)+entry
        all_spectra_2.append(tupl)

mzml_location, ids, mz, intensities = zip(*all_spectra_2)

spectra_df_2 = pd.DataFrame({'file':mzml_location,'id':ids,'mz':mz,'intensities':intensities})

#### MERGE: mzid + mzml

merged_df_2 = pd.merge(mzid_df_2,spectra_df_2,how='left',on=['file','id'])

merged_df_2 = merged_df_2[['id','seq','mz','intensities']]

#%%
hdf.put(value=merged_df_2, key="df2")

#%%
os.chdir('./train_3')

mzid_files_3 = glob.glob('*.mzid')
indexed_mzid_3 = mzid.chain.from_iterable(mzid_files_3, use_index=True)


def _parse_mzid_entry(entry):
    spectrum_id = str(entry['spectrumID'])
    seq = str(entry['SpectrumIdentificationItem'][0]['PeptideSequence'])
    try: 
        mods = entry['SpectrumIdentificationItem'][0]['Modification'] 
    except:
        mods = None
    rank = int(entry['SpectrumIdentificationItem'][0]['rank'])
    file_location = str(entry['name'])
    return file_location,spectrum_id,seq,mods,rank

all_mzid_3 = []
for entry in(indexed_mzid_3):
    all_mzid_3.append(_parse_mzid_entry(entry))
  

file_location,spectrum_ids,seq,mods,rank = zip(*all_mzid_3)

mzid_df_3 = pd.DataFrame({'file':file_location,'id':spectrum_ids,'seq':seq})

def _parse_mzml_entry(entry):
    ID = str(entry['id'])
    mz = np.array(entry['m/z array'])
    intensities = np.array(entry['intensity array'])
    return ID, mz, intensities

all_spectra_3 = []

for file in np.unique(file_location):
    print(file)
    indexed = mzml.MzML(file)
    for i,entry in enumerate(indexed.map(_parse_mzml_entry)):
        tupl = (file,)+entry
        all_spectra_3.append(tupl)

mzml_location, ids, mz, intensities = zip(*all_spectra_3)

spectra_df_3 = pd.DataFrame({'file':mzml_location,'id':ids,'mz':mz,'intensities':intensities})

#### MERGE: mzid + mzml

merged_df_3 = pd.merge(mzid_df_3,spectra_df_3,how='left',on=['file','id'])

merged_df_3 = merged_df_3[['id','seq','mz','intensities']]

#%%
hdf.put(value=merged_df_3, key="df3")

#%%
os.chdir('./train_4')

mzid_files_4 = glob.glob('*.mzid')
indexed_mzid_4 = mzid.chain.from_iterable(mzid_files_4, use_index=True)


def _parse_mzid_entry(entry):
    spectrum_id = str(entry['spectrumID'])
    seq = str(entry['SpectrumIdentificationItem'][0]['PeptideSequence'])
    try: 
        mods = entry['SpectrumIdentificationItem'][0]['Modification'] 
    except:
        mods = None
    rank = int(entry['SpectrumIdentificationItem'][0]['rank'])
    file_location = str(entry['name'])
    return file_location,spectrum_id,seq,mods,rank

all_mzid_4 = []
for entry in(indexed_mzid_4):
    all_mzid_4.append(_parse_mzid_entry(entry))
  

file_location,spectrum_ids,seq,mods,rank = zip(*all_mzid_4)

mzid_df_4 = pd.DataFrame({'file':file_location,'id':spectrum_ids,'seq':seq})

def _parse_mzml_entry(entry):
    ID = str(entry['id'])
    mz = np.array(entry['m/z array'])
    intensities = np.array(entry['intensity array'])
    return ID, mz, intensities

all_spectra_4 = []

for file in np.unique(file_location):
    print(file)
    indexed = mzml.MzML(file)
    for i,entry in enumerate(indexed.map(_parse_mzml_entry)):
        tupl = (file,)+entry
        all_spectra_4.append(tupl)

mzml_location, ids, mz, intensities = zip(*all_spectra_4)

spectra_df_4= pd.DataFrame({'file':mzml_location,'id':ids,'mz':mz,'intensities':intensities})

#### MERGE: mzid + mzml

merged_df_4 = pd.merge(mzid_df_4,spectra_df_4,how='left',on=['file','id'])

merged_df_4 = merged_df_4[['id','seq','mz','intensities']]


#%%
hdf.put(value=merged_df_4, key="df4")

#%%
os.chdir('./train_5')

mzid_files_5 = glob.glob('*.mzid')
indexed_mzid_5 = mzid.chain.from_iterable(mzid_files_5, use_index=True)


def _parse_mzid_entry(entry):
    spectrum_id = str(entry['spectrumID'])
    seq = str(entry['SpectrumIdentificationItem'][0]['PeptideSequence'])
    try: 
        mods = entry['SpectrumIdentificationItem'][0]['Modification'] 
    except:
        mods = None
    rank = int(entry['SpectrumIdentificationItem'][0]['rank'])
    file_location = str(entry['name'])
    return file_location,spectrum_id,seq,mods,rank

all_mzid_5 = []
for entry in(indexed_mzid_5):
    all_mzid_5.append(_parse_mzid_entry(entry))
  

file_location,spectrum_ids,seq,mods,rank = zip(*all_mzid_5)

mzid_df_5 = pd.DataFrame({'file':file_location,'id':spectrum_ids,'seq':seq})

def _parse_mzml_entry(entry):
    ID = str(entry['id'])
    mz = np.array(entry['m/z array'])
    intensities = np.array(entry['intensity array'])
    return ID, mz, intensities

all_spectra_5 = []

for file in np.unique(file_location):
    print(file)
    indexed = mzml.MzML(file)
    for i,entry in enumerate(indexed.map(_parse_mzml_entry)):
        tupl = (file,)+entry
        all_spectra_5.append(tupl)

mzml_location, ids, mz, intensities = zip(*all_spectra_5)

spectra_df_5 = pd.DataFrame({'file':mzml_location,'id':ids,'mz':mz,'intensities':intensities})

#### MERGE: mzid + mzml

merged_df_5 = pd.merge(mzid_df_5,spectra_df_5,how='left',on=['file','id'])

merged_df_5 = merged_df_5[['id','seq','mz','intensities']]

#%%
hdf.put(value=merged_df_5, key="df5")


#%%
os.chdir('./train_6')

mzid_files_6 = glob.glob('*.mzid')
indexed_mzid_6 = mzid.chain.from_iterable(mzid_files_6, use_index=True)


def _parse_mzid_entry(entry):
    spectrum_id = str(entry['spectrumID'])
    seq = str(entry['SpectrumIdentificationItem'][0]['PeptideSequence'])
    try: 
        mods = entry['SpectrumIdentificationItem'][0]['Modification'] 
    except:
        mods = None
    rank = int(entry['SpectrumIdentificationItem'][0]['rank'])
    file_location = str(entry['name'])
    return file_location,spectrum_id,seq,mods,rank

all_mzid_6 = []
for entry in(indexed_mzid_6):
    all_mzid_6.append(_parse_mzid_entry(entry))
  

file_location,spectrum_ids,seq,mods,rank = zip(*all_mzid_6)

mzid_df_6 = pd.DataFrame({'file':file_location,'id':spectrum_ids,'seq':seq})

def _parse_mzml_entry(entry):
    ID = str(entry['id'])
    mz = np.array(entry['m/z array'])
    intensities = np.array(entry['intensity array'])
    return ID, mz, intensities

all_spectra_6 = []

for file in np.unique(file_location):
    print(file)
    indexed = mzml.MzML(file)
    for i,entry in enumerate(indexed.map(_parse_mzml_entry)):
        tupl = (file,)+entry
        all_spectra_6.append(tupl)

mzml_location, ids, mz, intensities = zip(*all_spectra_6)

spectra_df_6 = pd.DataFrame({'file':mzml_location,'id':ids,'mz':mz,'intensities':intensities})

#### MERGE: mzid + mzml

merged_df_6 = pd.merge(mzid_df_6,spectra_df_6,how='left',on=['file','id'])

merged_df_6 = merged_df_6[['id','seq','mz','intensities']]

#%%
hdf.put(value=merged_df_6, key="df6")


#%%
os.chdir('./train_7')

mzid_files_7 = glob.glob('*.mzid')
indexed_mzid_7 = mzid.chain.from_iterable(mzid_files_7, use_index=True)


def _parse_mzid_entry(entry):
    spectrum_id = str(entry['spectrumID'])
    seq = str(entry['SpectrumIdentificationItem'][0]['PeptideSequence'])
    try: 
        mods = entry['SpectrumIdentificationItem'][0]['Modification'] 
    except:
        mods = None
    rank = int(entry['SpectrumIdentificationItem'][0]['rank'])
    file_location = str(entry['name'])
    return file_location,spectrum_id,seq,mods,rank

all_mzid_7 = []
for entry in(indexed_mzid_7):
    all_mzid_7.append(_parse_mzid_entry(entry))
  

file_location,spectrum_ids,seq,mods,rank = zip(*all_mzid_7)

mzid_df_7 = pd.DataFrame({'file':file_location,'id':spectrum_ids,'seq':seq})

def _parse_mzml_entry(entry):
    ID = str(entry['id'])
    mz = np.array(entry['m/z array'])
    intensities = np.array(entry['intensity array'])
    return ID, mz, intensities

all_spectra_7 = []

for file in np.unique(file_location):
    print(file)
    indexed = mzml.MzML(file)
    for i,entry in enumerate(indexed.map(_parse_mzml_entry)):
        tupl = (file,)+entry
        all_spectra_7.append(tupl)

mzml_location, ids, mz, intensities = zip(*all_spectra_7)

spectra_df_7 = pd.DataFrame({'file':mzml_location,'id':ids,'mz':mz,'intensities':intensities})


#### MERGE: mzid + mzml

merged_df_7 = pd.merge(mzid_df_7,spectra_df_7,how='left',on=['file','id'])

merged_df_7 = merged_df_7[['id','seq','mz','intensities']]


#%%
hdf.put(value=merged_df_7, key="df7")

#%%
os.chdir('./train_8')

mzid_files_8 = glob.glob('*.mzid')
indexed_mzid_8 = mzid.chain.from_iterable(mzid_files_8, use_index=True)


def _parse_mzid_entry(entry):
    spectrum_id = str(entry['spectrumID'])
    seq = str(entry['SpectrumIdentificationItem'][0]['PeptideSequence'])
    try: 
        mods = entry['SpectrumIdentificationItem'][0]['Modification'] 
    except:
        mods = None
    rank = int(entry['SpectrumIdentificationItem'][0]['rank'])
    file_location = str(entry['name'])
    return file_location,spectrum_id,seq,mods,rank

all_mzid_8 = []
for entry in(indexed_mzid_8):
    all_mzid_8.append(_parse_mzid_entry(entry))
  

file_location,spectrum_ids,seq,mods,rank = zip(*all_mzid_8)

mzid_df_8 = pd.DataFrame({'file':file_location,'id':spectrum_ids,'seq':seq})

def _parse_mzml_entry(entry):
    ID = str(entry['id'])
    mz = np.array(entry['m/z array'])
    intensities = np.array(entry['intensity array'])
    return ID, mz, intensities

all_spectra_8 = []

for file in np.unique(file_location):
    print(file)
    indexed = mzml.MzML(file)
    for i,entry in enumerate(indexed.map(_parse_mzml_entry)):
        tupl = (file,)+entry
        all_spectra_8.append(tupl)

mzml_location, ids, mz, intensities = zip(*all_spectra_8)

spectra_df_8 = pd.DataFrame({'file':mzml_location,'id':ids,'mz':mz,'intensities':intensities})


#### MERGE: mzid + mzml

merged_df_8 = pd.merge(mzid_df_8,spectra_df_8,how='left',on=['file','id'])

merged_df_8 = merged_df_8[['id','seq','mz','intensities']]

#%%
hdf.put(value=merged_df_8, key="df8")

#%%
os.chdir('./train_9')

mzid_files_9 = glob.glob('*.mzid')
indexed_mzid_9 = mzid.chain.from_iterable(mzid_files_9, use_index=True)


def _parse_mzid_entry(entry):
    spectrum_id = str(entry['spectrumID'])
    seq = str(entry['SpectrumIdentificationItem'][0]['PeptideSequence'])
    try: 
        mods = entry['SpectrumIdentificationItem'][0]['Modification'] 
    except:
        mods = None
    rank = int(entry['SpectrumIdentificationItem'][0]['rank'])
    file_location = str(entry['name'])
    return file_location,spectrum_id,seq,mods,rank

all_mzid_9 = []
for entry in(indexed_mzid_9):
    all_mzid_9.append(_parse_mzid_entry(entry))
  

file_location,spectrum_ids,seq,mods,rank = zip(*all_mzid_9)

mzid_df_9 = pd.DataFrame({'file':file_location,'id':spectrum_ids,'seq':seq})

def _parse_mzml_entry(entry):
    ID = str(entry['id'])
    mz = np.array(entry['m/z array'])
    intensities = np.array(entry['intensity array'])
    return ID, mz, intensities

all_spectra_9 = []

for file in np.unique(file_location):
    print(file)
    indexed = mzml.MzML(file)
    for i,entry in enumerate(indexed.map(_parse_mzml_entry)):
        tupl = (file,)+entry
        all_spectra_9.append(tupl)

mzml_location, ids, mz, intensities = zip(*all_spectra_9)

spectra_df_9 = pd.DataFrame({'file':mzml_location,'id':ids,'mz':mz,'intensities':intensities})


#### MERGE: mzid + mzml

merged_df_9 = pd.merge(mzid_df_9,spectra_df_9,how='left',on=['file','id'])

merged_df_9 = merged_df_9[['id','seq','mz','intensities']]

#%%
hdf.put(value=merged_df_9, key="df9")

#%%
os.chdir('./train_10')

mzid_files_10 = glob.glob('*.mzid')
indexed_mzid_10 = mzid.chain.from_iterable(mzid_files_10, use_index=True)


def _parse_mzid_entry(entry):
    spectrum_id = str(entry['spectrumID'])
    seq = str(entry['SpectrumIdentificationItem'][0]['PeptideSequence'])
    try: 
        mods = entry['SpectrumIdentificationItem'][0]['Modification'] 
    except:
        mods = None
    rank = int(entry['SpectrumIdentificationItem'][0]['rank'])
    file_location = str(entry['name'])
    return file_location,spectrum_id,seq,mods,rank

all_mzid_10 = []
for entry in(indexed_mzid_10):
    all_mzid_10.append(_parse_mzid_entry(entry))
  

file_location,spectrum_ids,seq,mods,rank = zip(*all_mzid_10)

mzid_df_10 = pd.DataFrame({'file':file_location,'id':spectrum_ids,'seq':seq})

def _parse_mzml_entry(entry):
    ID = str(entry['id'])
    mz = np.array(entry['m/z array'])
    intensities = np.array(entry['intensity array'])
    return ID, mz, intensities

all_spectra_10 = []

for file in np.unique(file_location):
    print(file)
    indexed = mzml.MzML(file)
    for i,entry in enumerate(indexed.map(_parse_mzml_entry)):
        tupl = (file,)+entry
        all_spectra_10.append(tupl)

mzml_location, ids, mz, intensities = zip(*all_spectra_10)

spectra_df_10 = pd.DataFrame({'file':mzml_location,'id':ids,'mz':mz,'intensities':intensities})


#### MERGE: mzid + mzml

merged_df_10 = pd.merge(mzid_df_10,spectra_df_10,how='left',on=['file','id'])

merged_df_10 = merged_df_10[['id','seq','mz','intensities']]

#%%
hdf.put(value=merged_df_10, key="df10")

#%%
os.chdir('./train_11')

mzid_files_11 = glob.glob('*.mzid')
indexed_mzid_11 = mzid.chain.from_iterable(mzid_files_11, use_index=True)


def _parse_mzid_entry(entry):
    spectrum_id = str(entry['spectrumID'])
    seq = str(entry['SpectrumIdentificationItem'][0]['PeptideSequence'])
    try: 
        mods = entry['SpectrumIdentificationItem'][0]['Modification'] 
    except:
        mods = None
    rank = int(entry['SpectrumIdentificationItem'][0]['rank'])
    file_location = str(entry['name'])
    return file_location,spectrum_id,seq,mods,rank

all_mzid_11 = []
for entry in(indexed_mzid_11):
    all_mzid_11.append(_parse_mzid_entry(entry))
  

file_location,spectrum_ids,seq,mods,rank = zip(*all_mzid_11)

mzid_df_11 = pd.DataFrame({'file':file_location,'id':spectrum_ids,'seq':seq})

def _parse_mzml_entry(entry):
    ID = str(entry['id'])
    mz = np.array(entry['m/z array'])
    intensities = np.array(entry['intensity array'])
    return ID, mz, intensities

all_spectra_11 = []

for file in np.unique(file_location):
    print(file)
    indexed = mzml.MzML(file)
    for i,entry in enumerate(indexed.map(_parse_mzml_entry)):
        tupl = (file,)+entry
        all_spectra_11.append(tupl)

mzml_location, ids, mz, intensities = zip(*all_spectra_11)

spectra_df_11 = pd.DataFrame({'file':mzml_location,'id':ids,'mz':mz,'intensities':intensities})


#### MERGE: mzid + mzml

merged_df_11 = pd.merge(mzid_df_11,spectra_df_11,how='left',on=['file','id'])

merged_df_11 = merged_df_11[['id','seq','mz','intensities']]


#%%
hdf.put(value=merged_df_11, key="df11")

#%%
os.chdir('./train_12')

mzid_files_12 = glob.glob('*.mzid')
indexed_mzid_12 = mzid.chain.from_iterable(mzid_files_12, use_index=True)


def _parse_mzid_entry(entry):
    spectrum_id = str(entry['spectrumID'])
    seq = str(entry['SpectrumIdentificationItem'][0]['PeptideSequence'])
    try: 
        mods = entry['SpectrumIdentificationItem'][0]['Modification'] 
    except:
        mods = None
    rank = int(entry['SpectrumIdentificationItem'][0]['rank'])
    file_location = str(entry['name'])
    return file_location,spectrum_id,seq,mods,rank

all_mzid_12 = []
for entry in(indexed_mzid_12):
    all_mzid_12.append(_parse_mzid_entry(entry))
  

file_location,spectrum_ids,seq,mods,rank = zip(*all_mzid_12)

mzid_df_12 = pd.DataFrame({'file':file_location,'id':spectrum_ids,'seq':seq})

def _parse_mzml_entry(entry):
    ID = str(entry['id'])
    mz = np.array(entry['m/z array'])
    intensities = np.array(entry['intensity array'])
    return ID, mz, intensities

all_spectra_12 = []

for file in np.unique(file_location):
    print(file)
    indexed = mzml.MzML(file)
    for i,entry in enumerate(indexed.map(_parse_mzml_entry)):
        tupl = (file,)+entry
        all_spectra_12.append(tupl)

mzml_location, ids, mz, intensities = zip(*all_spectra_12)

spectra_df_12 = pd.DataFrame({'file':mzml_location,'id':ids,'mz':mz,'intensities':intensities})

merged_df_12 = pd.merge(mzid_df_12,spectra_df_12,how='left',on=['file','id'])

merged_df_12 = merged_df_12[['id','seq','mz','intensities']]

# %%
os.chdir('./train_13')

mzid_files_13 = glob.glob('*.mzid')
indexed_mzid_13 = mzid.chain.from_iterable(mzid_files_13, use_index=True)


def _parse_mzid_entry(entry):
    spectrum_id = str(entry['spectrumID'])
    seq = str(entry['SpectrumIdentificationItem'][0]['PeptideSequence'])
    try: 
        mods = entry['SpectrumIdentificationItem'][0]['Modification'] 
    except:
        mods = None
    rank = int(entry['SpectrumIdentificationItem'][0]['rank'])
    file_location = str(entry['name'])
    return file_location,spectrum_id,seq,mods,rank

all_mzid_13 = []
for entry in(indexed_mzid_13):
    all_mzid_13.append(_parse_mzid_entry(entry))
  

file_location,spectrum_ids,seq,mods,rank = zip(*all_mzid_13)

mzid_df_13 = pd.DataFrame({'file':file_location,'id':spectrum_ids,'seq':seq})

def _parse_mzml_entry(entry):
    ID = str(entry['id'])
    mz = np.array(entry['m/z array'])
    intensities = np.array(entry['intensity array'])
    return ID, mz, intensities

all_spectra_13 = []

for file in np.unique(file_location):
    print(file)
    indexed = mzml.MzML(file)
    for i,entry in enumerate(indexed.map(_parse_mzml_entry)):
        tupl = (file,)+entry
        all_spectra_13.append(tupl)

mzml_location, ids, mz, intensities = zip(*all_spectra_13)

spectra_df_13 = pd.DataFrame({'file':mzml_location,'id':ids,'mz':mz,'intensities':intensities})

merged_df_13 = pd.merge(mzid_df_13,spectra_df_13,how='left',on=['file','id'])

merged_df_13 = merged_df_13[['id','seq','mz','intensities']]

#%%
hdf.put(value=merged_df_13, key="df13")


# %%
os.chdir('./train_14')

mzid_files_14 = glob.glob('*.mzid')
indexed_mzid_14 = mzid.chain.from_iterable(mzid_files_14, use_index=True)


def _parse_mzid_entry(entry):
    spectrum_id = str(entry['spectrumID'])
    seq = str(entry['SpectrumIdentificationItem'][0]['PeptideSequence'])
    try: 
        mods = entry['SpectrumIdentificationItem'][0]['Modification'] 
    except:
        mods = None
    rank = int(entry['SpectrumIdentificationItem'][0]['rank'])
    file_location = str(entry['name'])
    return file_location,spectrum_id,seq,mods,rank

all_mzid_14 = []
for entry in(indexed_mzid_14):
    all_mzid_14.append(_parse_mzid_entry(entry))
  

file_location,spectrum_ids,seq,mods,rank = zip(*all_mzid_14)

mzid_df_14 = pd.DataFrame({'file':file_location,'id':spectrum_ids,'seq':seq})

def _parse_mzml_entry(entry):
    ID = str(entry['id'])
    mz = np.array(entry['m/z array'])
    intensities = np.array(entry['intensity array'])
    return ID, mz, intensities

all_spectra_14 = []

for file in np.unique(file_location):
    print(file)
    indexed = mzml.MzML(file)
    for i,entry in enumerate(indexed.map(_parse_mzml_entry)):
        tupl = (file,)+entry
        all_spectra_14.append(tupl)

mzml_location, ids, mz, intensities = zip(*all_spectra_14)

spectra_df_14 = pd.DataFrame({'file':mzml_location,'id':ids,'mz':mz,'intensities':intensities})

merged_df_14 = pd.merge(mzid_df_14,spectra_df_14,how='left',on=['file','id'])

merged_df_14 = merged_df_14[['id','seq','mz','intensities']]

#%%
hdf.put(value=merged_df_14, key="df14")

    
# %%
os.chdir('./train_15')

mzid_files_15 = glob.glob('*.mzid')
indexed_mzid_15 = mzid.chain.from_iterable(mzid_files_15, use_index=True)


def _parse_mzid_entry(entry):
    spectrum_id = str(entry['spectrumID'])
    seq = str(entry['SpectrumIdentificationItem'][0]['PeptideSequence'])
    try: 
        mods = entry['SpectrumIdentificationItem'][0]['Modification'] 
    except:
        mods = None
    rank = int(entry['SpectrumIdentificationItem'][0]['rank'])
    file_location = str(entry['name'])
    return file_location,spectrum_id,seq,mods,rank

all_mzid_15 = []
for entry in(indexed_mzid_15):
    all_mzid_15.append(_parse_mzid_entry(entry))
  

file_location,spectrum_ids,seq,mods,rank = zip(*all_mzid_15)

mzid_df_15 = pd.DataFrame({'file':file_location,'id':spectrum_ids,'seq':seq})

def _parse_mzml_entry(entry):
    ID = str(entry['id'])
    mz = np.array(entry['m/z array'])
    intensities = np.array(entry['intensity array'])
    return ID, mz, intensities

all_spectra_15 = []

for file in np.unique(file_location):
    print(file)
    indexed = mzml.MzML(file)
    for i,entry in enumerate(indexed.map(_parse_mzml_entry)):
        tupl = (file,)+entry
        all_spectra_15.append(tupl)

mzml_location, ids, mz, intensities = zip(*all_spectra_15)

spectra_df_15 = pd.DataFrame({'file':mzml_location,'id':ids,'mz':mz,'intensities':intensities})

merged_df_15 = pd.merge(mzid_df_15,spectra_df_15,how='left',on=['file','id'])

merged_df_15 = merged_df_15[['id','seq','mz','intensities']]

#%%

hdf.put(value=merged_df_15, key="df15")

# %%
os.chdir('./train_16')

mzid_files_16 = glob.glob('*.mzid')
indexed_mzid_16 = mzid.chain.from_iterable(mzid_files_16, use_index=True)


def _parse_mzid_entry(entry):
    spectrum_id = str(entry['spectrumID'])
    seq = str(entry['SpectrumIdentificationItem'][0]['PeptideSequence'])
    try: 
        mods = entry['SpectrumIdentificationItem'][0]['Modification'] 
    except:
        mods = None
    rank = int(entry['SpectrumIdentificationItem'][0]['rank'])
    file_location = str(entry['name'])
    return file_location,spectrum_id,seq,mods,rank

all_mzid_16 = []
for entry in(indexed_mzid_16):
    all_mzid_16.append(_parse_mzid_entry(entry))
  

file_location,spectrum_ids,seq,mods,rank = zip(*all_mzid_16)

mzid_df_16 = pd.DataFrame({'file':file_location,'id':spectrum_ids,'seq':seq})

def _parse_mzml_entry(entry):
    ID = str(entry['id'])
    mz = np.array(entry['m/z array'])
    intensities = np.array(entry['intensity array'])
    return ID, mz, intensities

all_spectra_16 = []

for file in np.unique(file_location):
    print(file)
    indexed = mzml.MzML(file)
    for i,entry in enumerate(indexed.map(_parse_mzml_entry)):
        tupl = (file,)+entry
        all_spectra_16.append(tupl)

mzml_location, ids, mz, intensities = zip(*all_spectra_16)

spectra_df_16 = pd.DataFrame({'file':mzml_location,'id':ids,'mz':mz,'intensities':intensities})

merged_df_16 = pd.merge(mzid_df_16,spectra_df_16,how='left',on=['file','id'])

merged_df_16 = merged_df_16[['id','seq','mz','intensities']]


# %%

hdf.put(value=merged_df_16, key="df16")

# %%
os.chdir('./train_17')

mzid_files_17 = glob.glob('*.mzid')
indexed_mzid_17 = mzid.chain.from_iterable(mzid_files_17, use_index=True)


def _parse_mzid_entry(entry):
    spectrum_id = str(entry['spectrumID'])
    seq = str(entry['SpectrumIdentificationItem'][0]['PeptideSequence'])
    try: 
        mods = entry['SpectrumIdentificationItem'][0]['Modification'] 
    except:
        mods = None
    rank = int(entry['SpectrumIdentificationItem'][0]['rank'])
    file_location = str(entry['name'])
    return file_location,spectrum_id,seq,mods,rank

all_mzid_17 = []
for entry in(indexed_mzid_17):
    all_mzid_17.append(_parse_mzid_entry(entry))
  

file_location,spectrum_ids,seq,mods,rank = zip(*all_mzid)

mzid_df_17 = pd.DataFrame({'file':file_location,'id':spectrum_ids,'seq':seq})

def _parse_mzml_entry(entry):
    ID = str(entry['id'])
    mz = np.array(entry['m/z array'])
    intensities = np.array(entry['intensity array'])
    return ID, mz, intensities

all_spectra_17 = []

for file in np.unique(file_location):
    print(file)
    indexed = mzml.MzML(file)
    for i,entry in enumerate(indexed.map(_parse_mzml_entry)):
        tupl = (file,)+entry
        all_spectra_17.append(tupl)

mzml_location, ids, mz, intensities = zip(*all_spectra_17)

spectra_df_17 = pd.DataFrame({'file':mzml_location,'id':ids,'mz':mz,'intensities':intensities})

merged_df_17 = pd.merge(mzid_df_17,spectra_df_17,how='left',on=['file','id'])

merged_df_17 = merged_df_17[['id','seq','mz','intensities']]
# %%
hdf.put(value=merged_df_17, key="df17")

