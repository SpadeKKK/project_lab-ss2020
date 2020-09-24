#%%
### Import dependencies for downloading and reading a file from the web
import os
import urllib.request
from urllib.request import urlretrieve

# %%
### Download NCBI-Small .fasta file
if not os.path.isfile("ncbismall.fasta"):
    print ("Downloading the NCBI Small FASTAfile ...")
    urlretrieve('https://owncloud.hpi.de/s/fa0aV3lp4Mu8Upq/download', 'ncbismall.fasta')
    print ("Done!")

# %%
### Import dependencies for manipulation of protein sequences
import Bio
from Bio import SeqIO
import pyteomics
from pyteomics import fasta, parser
import pandas as pd


### Extract relevant biological information from the header of each protein refseq
protein_id = []
seq_length = []
proteinName = []
specieName = []

peptide_id = []
peptide_list= []

for i, seq_record in enumerate(SeqIO.parse("ncbismall.fasta", "fasta")):
    protein_id.append(seq_record.id)
    seq_length.append(len(seq_record))
    specieName.append(seq_record.description.split("[")[1])
    prot = seq_record.description.split("[")[0]
    protName = " ".join(prot.split(" ")[1:])
    proteinName.append(protName)


    ### Cleave protein sequence using trypsin enzyme, into strings of peptide sequences
    new_peptides = parser.cleave(str(seq_record.seq), parser.expasy_rules['trypsin'])
    for peptide in new_peptides:
        ### Select peptide sequences of lengths 7 to 24
        if len(peptide) > 24 or len(peptide) < 7:
            continue

        peptide_id.append(seq_record.id)
        peptide_list.append(peptide)
        
#%%
print(i)
### Number of proteins = 1699847
#%%
### Create a comprehensive dictionary containing the biologial information from each protein
ncbi_dict = {
    'pid': protein_id,
    'name_of_protein':proteinName,
    'name_of_specie': specieName,
    'length_of_peptide': seq_length
    }

### Convert the dictionary to a Dataframe
ncbi_df = pd.DataFrame(ncbi_dict)

### Save the Dataframe as a .csv file
ncbi_df.to_csv('ncbiDFsmall.csv')

#%%
### Create a Dataframe that maps each Protein to their Peptide sequence
ncbi_protpep = {'peptide_id': peptide_id,
                'peptide_list': peptide_list
                }

ncbi_protpep_df = pd.DataFrame(ncbi_protpep)
print("Number of cleaved peptides: ", len(ncbi_protpep_df))
### Number of cleaved peptides = 25627419

ncbi_protpep_df = ncbi_protpep_df[['peptide_list', 'peptide_id']]

### Convert dataframe to list
ncbi_protpep_list = list(ncbi_protpep_df.values.tolist())

### Visualize the list arrangement
#ncbi_protpep_list[0]        #prints the list/pair
#ncbi_protpep_list[0][0]     #prints peptide; [0][1] prints the ID


#%%
### Create a Dictionary of unique peptides and mapped to their protein IDs
ncbi_peptide_dict = {}
for a in ncbi_protpep_list:
    if a[0] in ncbi_peptide_dict:
        s = ncbi_peptide_dict[a[0]]
        s.append(a[1])
        ncbi_peptide_dict[a[0]] = s
    else:
        ncbi_peptide_dict[a[0]] = [a[1]]

print(len(ncbi_peptide_dict.keys()))
#ncbi_peptide_dict

# %%
### Save unique peptide-protein as dictionary
import csv

with open('ncbi_uniquepeptides.csv', 'w') as f:
    for key in ncbi_peptide_dict.keys():
        f.write("%s, %s\n" % (key, ncbi_peptide_dict[key]))


# %%
### Filter the dataframe, get rid of X, B, J, U & Z alphabets
### Overwrite with the filtered dataset Unique peptide dictionary

ncbi_protpep_df_X = ncbi_protpep_df[~ncbi_protpep_df.peptide_list.str.contains("X")]

ncbi_protpep_df_B = ncbi_protpep_df_X[~ncbi_protpep_df_X.peptide_list.str.contains("B")]

ncbi_protpep_df_J = ncbi_protpep_df_B[~ncbi_protpep_df_B.peptide_list.str.contains("J")]

ncbi_protpep_df_U = ncbi_protpep_df_J[~ncbi_protpep_df_J.peptide_list.str.contains("U")]

ncbi_protpep_df_Z = ncbi_protpep_df_U[~ncbi_protpep_df_U.peptide_list.str.contains("Z")]

#%%
### Overwrite with the new list of peptide-protein pairs

ncbismall_list = list(ncbi_protpep_df_Z.values.tolist())

# %%
### Recreate dictionary with unique peptides
ncbi_unique_dict = {}
for a in ncbismall_list:
    if a[0] in ncbi_unique_dict:
        s = ncbi_unique_dict[a[0]]
        s.append(a[1])
        ncbi_unique_dict[a[0]] = s
    else:
        ncbi_unique_dict[a[0]] = [a[1]]

print(len(ncbi_unique_dict.keys()))
ncbi_unique_dict

#%%
### Print final length of unique peptides
unique_peptides = set(ncbi_protpep_df_Z.loc[:, 'peptide_list'])
print('Final Length of Unique Peptides: ', len(unique_peptides))
### 6780082

# %%
### Create a dataframe which maps unique peptides to multiple protein ids
unique_dt = pd.DataFrame(ncbi_unique_dict.items(), columns=['peptide', 'protein_id'])
unique_dt
# %%
### Unique peptides is saved into a file;it was better to do so because of its huge size when reading.
# The file can be read from the path provided as below: 
# unique_dt = pd.read_hdf('/home/ubuntu/data/jiahao/unique_df.hdf5', key="df", mode="r")