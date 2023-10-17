#!/usr/bin/env python
# coding: utf-8

# # Supplementary Notebook 2: Correlation between PTR and abundance
# ## Paper: Novel Approach for Microbiome Analysis Using Bacterial Replication Rates and Causal Inference to Determine Resistome Potential
# ### Vitalii Stebliankin, Musfiqur Sazal, Camilo Valdes, Kalai Mathee, and GiriNarasimhan
# 
# #### Dataset: Gibson et al. (BioProject ID: PRJNA301903)
# 
# ### Get correlation between PTR and abundance
# (reffered to Fig. 2A in the main manuscript)

# In[1]:


import pandas as pd
from statsmodels.stats.multitest import multipletests
from scipy.stats import spearmanr
import os

#intermediate_dir = "{}/intermediate_files".format(out_dir)
#if not os.path.exists(out_dir):
#    os.mkdir(out_dir)
# if not os.path.exists(intermediate_dir):
#     os.mkdir(intermediate_dir)



def get_species_list(ptr_df):
    columns = ptr_df.columns
    species_list=[]
    for col in columns:
        if "PTR" in col:
            species = col.replace("#PTR", "")
            species_list.append(species)
    print("Total of {} species in the dataset.".format(len(species_list)))
    return species_list

class CorrelationPTRPlugin:
 def input(self, infile):
    self.PTR_file = infile#"analysis-out/1-FilteringPTR/PTR_species_filtered_metadata_major.csv"

 def run(self):
     pass

 def output(self, outfile):
     out_dir = outfile#"analysis-out/2-CorrelationPTR-Abundance"
     ptr_df = pd.read_csv(self.PTR_file, index_col=0)
     out_file = out_dir+"/PTR_abundance_corr.csv"
     species_list = get_species_list(ptr_df)
     species_list_PTR = [x+"#PTR" for x in species_list]
     # ptr_only_df = ptr_df[species_list_PTR]
     # ptr_only_df.columns = ptr_only_df.columns.str.replace("#PTR","")

     ptr_dict = {"species":[],"spearmanr":[], "pval":[], "averageAbundance":[]}
     for species in species_list:
      # Drop NA values:
      tmp_df = ptr_df[ptr_df[species+"#PTR"].notnull()]
      ptr_dict["species"].append(species)
      if len(tmp_df)>0:
        average_abundance = ptr_df[species+"#abundance"].mean()
        r, pval = spearmanr(tmp_df[species+"#PTR"], tmp_df[species+"#abundance"])
      else:
        r, pval = 0, 1
        average_abundance = 0
      ptr_dict["spearmanr"].append(r)
      ptr_dict["pval"].append(pval)
      ptr_dict["averageAbundance"].append(average_abundance)
     df = pd.DataFrame(ptr_dict)
     df = df[["species", "spearmanr", "pval", "averageAbundance"]]
     df.to_csv(out_file)
     # Get average abundance of every species
     print(df)


