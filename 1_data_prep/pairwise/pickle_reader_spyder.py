import pandas as pd
import dill
import numpy as np
import seaborn as sns

test = np.load("/media/sf_Shared_Folder/PhD_work/data/scz/UCLA/pydata//AROMA_2P/sub-10206.npy")
pkl_file="/media/sf_Shared_Folder/PhD_work/data/scz/UCLA/pydata//AROMA_2P/sub-10206/calc.pkl"


with open(pkl_file,'rb') as f:
  calc = dill.load(f)

SPI="pli_multitaper_max_fs-1_fmin-0_fmax-0-5"
res = calc.table[SPI]
ax = sns.heatmap(res, annot=False)
ax.invert_yaxis()
  
sample_dict = {}

  
# Iterate over each SPI
for SPI in calc._spis:
    res = calc.table[SPI]
    res.reset_index(level=0, inplace=True)
    res = res.rename(columns={"index": "brain_region_1"})
    
    # Convert dataframe from wide to long
    res_long = pd.melt(res, id_vars='brain_region_1')
    res_long = res_long.rename(columns={"process": "brain_region_2"})
    
    # Set SPI and subject ID
    res_long["SPI"] = SPI
    
    res_long["brain_region_1"] = res_long["brain_region_1"].replace("proc-", "", regex=True).astype(int) + 1
    res_long["brain_region_2"] = res_long["brain_region_2"].replace("proc-", "", regex=True).astype(int) + 1
    
    # Omit pairs where the two brain regions are the same
    res_long = res_long[res_long["brain_region_1"] != res_long["brain_region_2"]]
    
    # add results to the subject's dictionary
    sample_dict[SPI] = res_long

# Combine dictionary of SPI dataframes into one dataframe
sample_df = pd.concat(sample_dict, axis=0)
