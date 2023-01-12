import os
import dill
import dplython as dpl
from string import capitalize
from dplython import (DplyFrame, 
                      X, 
                      diamonds, 
                      dfilter,
                      select,
                      rename,
                      sift, 
                      sample_n,
                      sample_frac, 
                      head, 
                      arrange,
                      mutate,
                      nrow,
                      group_by,
                      summarize, 
                      DelayFunction) 
import pandas as pd
import numpy as np

# Define paths specific to the UCLA CNP and ABIDE ASD datasets
univariate_feature_set = "catch22"
pairwise_feature_set = "pyspi14"
github_dir = "/Users/abry4213/github/fMRI_FeaturesDisorders/"
data_path = "/Users/abry4213//data/UCLA_CNP_ABIDE_ASD/"
rdata_path = data_path + "processed_data/Rdata/"
pydata_path = data_path + "processed_data/pydata/"
UCLA_CNP_noise_proc = "AROMA+2P+GMR"
ABIDE_ASD_noise_proc = "FC1000"


############################## Metadata ########################################

if not os.exists(data_path + "study_metadata/UCLA_CNP_ABIDE_ASD_sample_metadata.pkl"):
    # Metadata for UCLA CNP
    UCLA_CNP_metadata = DplyFrame(pd.read_csv(data_path + "study_metadata/UCLA_CNP_participants.csv"))
    UCLA_CNP_metadata = (UCLA_CNP_metadata >> 
                          select(X.participant_id, X.diagnosis, X.age, X.gender) >> 
                          mutate(Sample_ID = X.participant_id,
                                 Diagnosis = X.diagnosis,
                                 Age = X.age,
                                 Sex = X.gender) >>
                          select(X.Sample_ID, X.Diagnosis, X.Age, X.Sex) >>
                          mutate(Study = "UCLA_CNP"))
    # Cast to pandas
    UCLA_CNP_metadata = pd.DataFrame(UCLA_CNP_metadata)
    
    # Metadata for ABIDE ASD
    lookup_table = ["Control", "ASD"]
    ABIDE_ASD_metadata = DplyFrame(pd.read_csv(data_path + "study_metadata/ABIDE_ASD_participants.csv"))
    ABIDE_ASD_metadata = (ABIDE_ASD_metadata >> 
                          select(X.subject_id, X.asd, X.age, X.sex, X.site) >> 
                          mutate(Sample_ID = X.subject_id,
                                 ASD = X.asd,
                                 Age = X.age,
                                 Sex = X.sex,
                                 Site = X.site) >>
                          select(X.Sample_ID, X.ASD, X.Age, X.Sex, X.Site) >>
                          mutate(Study = "ABIDE_ASD"))
    
    ABIDE_ASD_metadata = pd.DataFrame(ABIDE_ASD_metadata)
    ABIDE_ASD_metadata["Diagnosis"] = [lookup_table[i] for i in ABIDE_ASD_metadata["ASD"].tolist()]
    ABIDE_ASD_metadata = ABIDE_ASD_metadata.drop(["ASD"], axis=1)
    
    # Merge all the metadata
    merged_metadata = pd.concat([UCLA_CNP_metadata, ABIDE_ASD_metadata])
    
    # Save to a pkl file
    with open(data_path + "study_metadata/UCLA_CNP_ABIDE_ASD_sample_metadata.pkl", "wb") as f:
        dill.dump(merged_metadata, f)

############################## Univariate ########################################

# catch22
if not os.exists(pydata_path + "study_metadata/UCLA_CNP_ABIDE_ASD_catch22.pkl"):
    UCLA_CNP_catch22 = 