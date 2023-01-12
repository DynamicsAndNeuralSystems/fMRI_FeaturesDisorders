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
data_path = "/Users/abry4213//data/UCLA_CNP/"
rdata_path = data_path + "processed_data/Rdata/"
pydata_path = data_path + "processed_data/pydata/"
noise_proc = "AROMA+2P+GMR"
metadata_CSV = data_path + "study_metadata/UCLA_CNP_participants.csv"

############################## Save to pickle ########################################

if not os.path.exists(data_path + "study_metadata/UCLA_CNP_sample_metadata.pkl"):
    # Metadata for UCLA CNP
    UCLA_CNP_metadata = DplyFrame(pd.read_csv(metadata_CSV))
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
    
    # Save to a pkl file
    with open(data_path + "study_metadata/UCLA_CNP_sample_metadata.pkl", "wb") as f:
        dill.dump(UCLA_CNP_metadata, f)
