#%%
import os
from pyteomics import mzid, mzml
import pandas as pd
import numpy as np
import glob

"""
Identically as how we did with the training data set, we randomly divided the test files into different
folders, then we generated different data frames and stored all of them in one single hdf file as our
validation daata set
"""

#%%
os.chdir('./test')

mzid_files=glob.glob('*.mzid')
indexed_mzid = mzid.chain.from_iterable(mzid_files,use_index=True)

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
hdf_test = pd.HDFStore('/home/ubuntu/data/jiahao/files/test.hdf5', mode='w')

#%%
hdf_test.put(value=merged_df, key="df")

#%%
os.chdir('./test_1')

mzid_files_1 = glob.glob('*.mzid')
indexed_mzid_1 = mzid.chain.from_iterable(mzid_files_1,use_index=True)

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
hdf_test.put(value=merged_df_1, key="df1")

#%%
os.chdir('./test_2')

mzid_files_2 = glob.glob('*.mzid')
indexed_mzid_2 = mzid.chain.from_iterable(mzid_files_2,use_index=True)

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

hdf_test.put(value=merged_df_2, key="df2")


#%%
os.chdir('./test_4')

mzid_files_4 = glob.glob('*.mzid')
indexed_mzid_4 = mzid.chain.from_iterable(mzid_files_4,use_index=True)

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

spectra_df_4 = pd.DataFrame({'file':mzml_location,'id':ids,'mz':mz,'intensities':intensities})

#### MERGE: mzid + mzml

merged_df_4 = pd.merge(mzid_df_4,spectra_df_4,how='left',on=['file','id'])

merged_df_4 = merged_df_4[['id','seq','mz','intensities']]

#%%

hdf_test.put(value=merged_df_4, key="df4")

#%%
os.chdir('./test_5')

mzid_files_5 = glob.glob('*.mzid')
indexed_mzid_5 = mzid.chain.from_iterable(mzid_files_5,use_index=True)

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

hdf_test.put(value=merged_df_5, key="df5")

#%%
os.chdir('./test_6')

mzid_files_6 = glob.glob('*.mzid')
indexed_mzid_6 = mzid.chain.from_iterable(mzid_files_6,use_index=True)

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
hdf_test.put(value=merged_df_6, key="df6")

# %%
os.chdir('./test_7')

mzid_files_7 = glob.glob('*.mzid')
indexed_mzid_7 = mzid.chain.from_iterable(mzid_files_7,use_index=True)

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
hdf_test.put(value=merged_df_7, key="df7")

# %%
os.chdir('./test_8')

mzid_files_8 = glob.glob('*.mzid')
indexed_mzid_8 = mzid.chain.from_iterable(mzid_files_8,use_index=True)

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
        all_spectra_7.append(tupl)

mzml_location, ids, mz, intensities = zip(*all_spectra_8)

spectra_df_8 = pd.DataFrame({'file':mzml_location,'id':ids,'mz':mz,'intensities':intensities})


merged_df_8 = pd.merge(mzid_df_8,spectra_df_8,how='left',on=['file','id'])

merged_df_8 = merged_df_8[['id','seq','mz','intensities']]


# %%

hdf.put(value=merged_df_8, key="df8")

