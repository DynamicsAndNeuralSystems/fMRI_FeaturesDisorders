#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 29 08:46:42 2023

@author: abry4213
"""


import os
import numpy as np
import nibabel as nib
import pandas as pd
from nilearn import plotting
from nilearn.image import math_img, resample_img, index_img, threshold_img
from matplotlib import pyplot as plt

# Define data path
input_data_path = "/Users/abry4213/data/UCLA_CNP/raw_data/neuroimaging/aparc_aseg_in_BOLD_space/"

# Helper function to parse 3D array
def process_array_data(neuroimg_array_3d, sample_ID):
    array_data_flattened = neuroimg_array_3d.flatten()
    
        
    # Get the indices of the original array dimensions
    indices = np.indices(neuroimg_array_3d.shape).reshape(neuroimg_array_3d.ndim, -1).T
    
    # Create a pandas DataFrame and filter to nonzero voxels only
    array_df = pd.DataFrame({
        'Value': array_data_flattened,
        'Voxel_x': indices[:, 0],
        'Voxel_y': indices[:, 1],
        'Voxel_z': indices[:, 2]
    }).query("Value != 0").reset_index()
    
    # Count voxels per ROI
    ROI_counts = (array_df
                  .groupby('Value').size()
                  .reset_index()
                  .rename(columns= {"Value": "ROI_Index",
                                    0: "Num_Voxels"})
                  .assign(Sample_ID = sample_ID)
                  )

    
    return(ROI_counts)

# Initialize a list to store ROI voxel counts
ROI_voxel_counts_list = []

# Iterate over nifti volumes
for aparc_aseg in os.listdir(input_data_path):
    sample_ID = aparc_aseg.replace("_bold_space-MNI152NLin2009cAsym_preproc_aparcaseg_roi.nii.gz", "")
    
    # Read in NIFTI volume data
    nifti_img = nib.load(f"{input_data_path}/{aparc_aseg}").get_fdata()
    
    # Calculate the number of voxels in each ROI
    ROI_voxel_counts = process_array_data(nifti_img, sample_ID)
    
    # Append the results to the growing list
    ROI_voxel_counts_list.append(ROI_voxel_counts)
    
# Concatenate results into one dataframe
ROI_voxel_counts = pd.concat(ROI_voxel_counts_list, axis=0).reset_index(drop=True)

# Same voxel volumes to a feather file
ROI_voxel_counts.to_feather("/Users/abry4213/data/UCLA_CNP/processed_data/aparc_aseg_BOLD_space_voxel_volumes.feather")
