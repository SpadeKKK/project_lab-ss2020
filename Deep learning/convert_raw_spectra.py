#%%
# import packages
from pyteomics import mzml
import numpy as np
import pandas as pd


# %%
# define parsing function 
def _parse_mzml_entry(entry):
    ID = str(entry['id'])
    mz = np.array(entry['m/z array'])
    intensities = np.array(entry['intensity array'])
    return ID, mz, intensities

all_spectra = []
data = '/home/ubuntu/data/jiahao/trp/output/Run1_U4_2000ng.mzML'
file = mzml.MzML(data)
for i,entry in enumerate(file.map(_parse_mzml_entry)):
    tupl = (data,)+entry
    all_spectra.append(tupl)

# %%
# generate pandas dataframe
mzml_location, ids, mz, intensities = zip(*all_spectra)

spectra_df = pd.DataFrame({'file':mzml_location,'id':ids,'mz':mz,'intensities':intensities})


# %%
#split the column id into 3
spectra_df[['controllertype', 'controllernumber', 'scan']] = spectra_df.id.str.split(expand=True)

# %%
#split scan into 2 columns to keep scan number in a single column
spectra_df[['SCAN', 'scan_number']] = spectra_df.scan.str.split("=", expand=True)

#%%
# save data frame
spectra_df = spectra_df[['file', 'id', 'mz', 'intensities', 'scan:number']]

spectra_df.to_hdf('spectra_for_knn_search.hdf5', key="df", mode="w")
