#!/usr/local/bin/python3

# Import required modules
import glob
import scipy.io as sio
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('ggplot')
import seaborn as sns
import sklearn
import analysisFunctions as af

from scipy.stats import zscore
from scipy import stats

# Read in and store the framewise displacement (fd) for the given dataset in a variable called fdAvgs,
# and create the TS_path_names and indices2Keep variables

# Store the fdAvgs and set a threshold fd
filePath = '/Users/AV/Dropbox/COBRE/movementData/fdAvgs_COBRE.txt'
fdAvgs = pd.read_csv(filePath,header=None);
threshold_fd = 0.5

# Store the path of the folder containing the subject data for the given dataset
subPath = '/Users/AV/Dropbox/COBRE/cfgData/'

# Need to alphabetise and store the subject file names into a variable
TS_path_names = sorted(glob.glob(subPath + '*.mat'))

# Filter the subjects based on their fd, and retain the subjects that have an fd < threshold_fd
TS_path_names, indices2Keep = af.removePathNames(filePath, threshold_fd, TS_path_names)
indices2Keep = indices2Keep.tolist()

# Adding 1 to every element in the array to convert to MATLAB indexing
indices2KeepMat = list(np.asarray(indices2Keep) + 1)

# Add a multi-level index to the tsData and store some key variables

element = 'element1_COBRE.txt' # Read in the feature matrix data from the saved .txt file
PyFeatList = 'PythonFeatureList.txt' # This text file contains the 22 feature names

# Add a multi-level index to the feature matrix and save into the variable, tsData
# Also store the number of ROIs and subjects in the data
tsData, ROIs, subjects, feats, featList = af.addIndices(element,subPath,PyFeatList)

# Choose which feature to analyse
feature = 1
featureName = featList[feature-1]

# Choose which ROI to analyse
ROI = 1

# Create the target column - unique for each dataset

# Select which dataset is being used
dataset = 'COBRE'

if dataset == 'UCLA':

    # Creating the target column
    targetCol = af.getTargetCol(TS_path_names)

elif dataset == 'COBRE':

    # Creating the target column
    csvPath = '/Users/AV/Dropbox/COBRE/participants.csv'
    COBRE = pd.read_csv(csvPath,header=None);

    targetCol = COBRE.iloc[1:,2]
    targetCol = targetCol.tolist()
    targetCol = pd.DataFrame(data=targetCol, columns=['target'])

    targetCol = targetCol.iloc[indices2Keep,:]
    targetCol = np.asarray(targetCol,dtype=np.int)

    # A '0' indicates a control subject and a '1' indicates a subject with SCZ
    targetColModified = np.where(targetCol==1, 0, targetCol) # First change the pre-existing 1s to 0s
    targetCol = np.where(targetCol==2, 1, targetColModified) # Then change the 2s to 1s

# Store and print the subject numbers within the dataset
Control, SCZ, Total, SCZ2Ctrl = af.giveMeSubjectNums(targetCol)

# Initialise a few boolean variables which decide what the outputs are
dispFigs = True
regAccOnly = False
featAccOnly = False

# Region by Region Analysis
af.Reg_by_Reg_Anal(ROI, tsData, targetCol, ROIs, indices2KeepMat, regAccOnly, dispFigs)

# Feature by Feature Analysis
af.Feat_by_Feat_Anal(feature, featureName, element, subPath, PyFeatList,
                     indices2KeepMat, targetCol, featAccOnly, dispFigs)

af.showMeTValHistograms(element,subPath,PyFeatList,dataset,indices2KeepMat,targetCol)

fp1 = '/Users/AV/Desktop/FeatureMatrixData/sortedFeats_UCLA.txt'
fp2 = '/Users/AV/Desktop/FeatureMatrixData/sortedFeats_COBRE.txt'

af.showMeJointPlot(fp1,fp2)

fp3 = '/Users/AV/Desktop/FeatureMatrixData/featBalancedAcc_UCLA_shuff.txt'
fp4 = '/Users/AV/Desktop/FeatureMatrixData/featBalancedAcc_COBRE_shuff.txt'

af.computeCombPVals(fp1,fp2,fp3,fp4,feats)

# Plot the fd vs balanced accuracy

# Import and store the fdArray
filePathA = '/Users/AV/Desktop/FeatureMatrixData/fdArray_COBRE_E1.txt'

af.giveMeFDvBalancedAcc(filePathA)

# Navigate to the 'fdAvgs_DATASET.txt' file
filePathB = '/Users/AV/Dropbox/COBRE/movementData/fdAvgs_COBRE.txt'

af.showMeFDAcrossSubs(filePathA,filePathB)
