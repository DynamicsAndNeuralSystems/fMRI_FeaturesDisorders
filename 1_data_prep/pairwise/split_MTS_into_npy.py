#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Load libraries
import numpy as np
import pandas as pd
import scipy.io
from scipy import stats
import csv

study = "/media/sf_Shared_Folder/PhD_work/"

###################### UCLA ##########################
# Load subject data
ucla_subj_info = pd.read_csv(study + "data/scz/UCLA/participants.csv", sep=",")

# Load data
mat_file = study + "data/scz/UCLA/new/UCLA_time_series_four_groups.mat"
ucla_data = scipy.io.loadmat(mat_file)

# Load preprocessing data
def flatten(t):
    return [item for sublist in t for item in sublist]
noise_proc_info = flatten([x.tolist() for x in ucla_data["noiseOptions"].tolist()[0]])

# ROI info
ROI_info = ucla_data["StructNames"].tolist()
ROI_info = [ROI.replace(" ", "") for ROI in ROI_info]
ROI_df = pd.DataFrame(data={"ROI": ROI_info,
                            "Index": list(range(1, len(ROI_info) + 1))})

# Write ROI index + name to a CSV
df.to_csv(study + "data/scz/UCLA/pydata/ROI_info.csv", sep=',',index=False)

# Load subject IDs
subject_IDs = flatten(flatten(ucla_data["subject_list"].tolist()))

# Read in time-series data and split by noise processing method (4th dimension)
TS_full = ucla_data["time_series"]
TS_data = np.split(TS_full, 3, axis=3)

# AROMA 2P
AROMA_2P_raw = TS_data[0][:,:,:,0]
AROMA_2P_split = np.split(AROMA_2P_raw, 260, axis=2)

def array_to_df_AROMA2P(i):
    data_array = AROMA_2P_split[i]
    array_2d = data_array[:,:,0]
    subjID = subject_IDs[i]
    num_times, num_regions = array_2d.shape
    print(f"Subject {subjID}: {num_regions} regions, {num_times} time points\n")
    data_norm = np.apply_along_axis(stats.zscore, 0, array_2d)
    data_norm = np.transpose(data_norm)
    np.save(study+f"data/scz/UCLA/pydata/AROMA_2P/{subjID}.npy", data_norm)

[array_to_df_AROMA2P(i) for i in list(range(len(AROMA_2P_split)))]

# AROMA 2P_GMR
AROMA_2P_GMR_raw = TS_data[1][:,:,:,0]
AROMA_2P_GMR_split = np.split(AROMA_2P_GMR_raw, 260, axis=2)

def array_to_df_AROMA2P_GMR(i):
    data_array = AROMA_2P_GMR_split[i]
    array_2d = data_array[:,:,0]
    subjID = subject_IDs[i]
    num_times, num_regions = array_2d.shape
    print(f"Subject {subjID}: {num_regions} regions, {num_times} time points\n")
    data_norm = np.apply_along_axis(stats.zscore, 0, array_2d)
    data_norm = np.transpose(data_norm)
    np.save(study+f"data/scz/UCLA/pydata/AROMA_2P_GMR/{subjID}.npy", data_norm)

[array_to_df_AROMA2P_GMR(i) for i in list(range(len(AROMA_2P_GMR_split)))]

# AROMA 2P_DiCER
AROMA_2P_DiCER_raw = TS_data[2][:,:,:,0]
AROMA_2P_DiCER_split = np.split(AROMA_2P_DiCER_raw, 260, axis=2)

def array_to_df_AROMA2P_DiCER(i):
    data_array = AROMA_2P_DiCER_split[i]
    array_2d = data_array[:,:,0]
    subjID = subject_IDs[i]
    num_times, num_regions = array_2d.shape
    print(f"Subject {subjID}: {num_regions} regions, {num_times} time points\n")
    data_norm = np.apply_along_axis(stats.zscore, 0, array_2d)
    data_norm = np.transpose(data_norm)
    np.save(study+f"data/scz/UCLA/pydata/AROMA_2P_DiCER/{subjID}.npy", data_norm)

[array_to_df_AROMA2P_DiCER(i) for i in list(range(len(AROMA_2P_DiCER_split)))]
