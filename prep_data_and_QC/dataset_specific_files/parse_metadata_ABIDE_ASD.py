import os.path
import dill
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

# Define paths specific to the UCLA CNP and ABIDE ASD datasets
univariate_feature_set = "catch22"
pairwise_feature_set = "pyspi14"
github_dir = "/Users/abry4213/github/fMRI_FeaturesDisorders/"
data_path = "/Users/abry4213//data/ABIDE_ASD/"
rdata_path = data_path + "processed_data/Rdata/"
pydata_path = data_path + "processed_data/pydata/"
noise_proc = "FC1000"

############################## Save to pickle ########################################

if not os.path.exists(data_path + "study_metadata/ABIDE_ASD_sample_metadata.pkl"):
    
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
    
    # Save to a pkl file
    with open(data_path + "study_metadata/ABIDE_ASD_sample_metadata.pkl", "wb") as f:
        dill.dump(ABIDE_ASD_metadata, f)
