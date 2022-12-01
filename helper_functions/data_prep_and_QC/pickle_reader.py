import pandas as pd
import dill

def extract_df_from_pkl(pkl_file):

    with open(pkl_file,'rb') as f:
      SPI_res = dill.load(f)

    # Iterate over each SPI
    SPI_res.columns = SPI_res.columns.to_flat_index()
    
    # Convert index to column
    SPI_res.reset_index(level=0, inplace=True)
    
     # Rename index as first brain region
    SPI_res = SPI_res.rename(columns={"index": "brain_region_1"})

    # Pivot data from wide to long
    SPI_res_long = pd.melt(SPI_res, id_vars="brain_region_1")
    SPI_res_long['SPI'], SPI_res_long['brain_region_2'] = SPI_res_long.variable.str
    
     # Remove variable column
    SPI_res_long = SPI_res_long.drop("variable", 1)
    
    # Return final dictionary
    return(SPI_res_long)

def read_in_df(pkl_file):

    with open(pkl_file,'rb') as f:
      SPI_res = dill.load(f)
    
    # Return dataframe
    return(SPI_res)