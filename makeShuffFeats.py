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
#-------------------------------------------------------------------------------

# Read in and store the framewise displacement (fd) for the given dataset in a variable called fdAvgs,
# and create the TS_path_names and indices2Keep variables

# Store the fdAvgs and set a threshold fd
filePath = '/Users/AV/Dropbox/COBRE/movementData/fdAvgs_COBRE.txt'
fdAvgs = pd.read_csv(filePath,header=None);
threshold_fd = 0.5


df3 = pd.DataFrame({'feature': [], 'featBalancedAcc': [], 'stdDev': [], 'svmWeights' : []})


# Store the path of the folder containing the subject data for the given dataset
subPath = '/Users/AV/Dropbox/COBRE/cfgData/'

# Need to alphabetise and store the subject file names into a variable
TS_path_names = sorted(glob.glob(subPath + '*.mat'))

# Filter the subjects based on their fd, and retain the subjects that have an fd < threshold_fd
TS_path_names, indices2Keep = af.removePathNames(filePath, threshold_fd, TS_path_names)
indices2Keep = indices2Keep.tolist()

# Adding 1 to every element in the array to convert to MATLAB indexing
indices2KeepMat = list(np.asarray(indices2Keep) + 1)

# print(indices2KeepMat)
#-------------------------------------------------------------------------------

# Add a multi-level index to the tsData and store some key variables

element = 'element1_COBRE.txt' # Read in the feature matrix data from the saved .txt file
PyFeatList = 'PythonFeatureList.txt' # This text file contains the 22 feature names

# Add a multi-level index to the feature matrix and save into the variable, tsData
# Also store the number of ROIs and subjects in the data
tsData, ROIs, subjects, feats, featList = af.addIndices(element,subPath,PyFeatList)
#-------------------------------------------------------------------------------

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
#-------------------------------------------------------------------------------

for feature in range(1, 23):

    # Choose which feature to analyse
    featureName = featList[feature-1]
    featSlice = af.getFeatSlice(ROIs,subjects,tsData,featureName,indices2KeepMat)
    DataSlice = featSlice
    DataSlice_zscored = DataSlice.apply(zscore)
    X = DataSlice_zscored
    y = np.ravel(targetCol)

    # Perform 10-Fold Cross Validation

    # Store the function's output as a variable
    scores = af.get10FoldCVScore(X,y)
    weights = giveMeSVMWeights(X,y)

    featBalancedAcc = scores.mean()
    stdDev = scores.std()

    df3 = df3.append({'feature': feature, 'featBalancedAcc': featBalancedAcc, 'stdDev': stdDev, 'svmWeights' : weights}, ignore_index=True)
    df3_sorted = df3.sort_values(by='featBalancedAcc',ascending=False)
    df3_sorted = df3_sorted.set_index('feature')

#     for i in range(1,1001):

#         iteration = i

#         # Shuffled target column
#         np.random.shuffle(targetCol)
#         y = np.ravel(targetCol)
#         #-------------------------------------------------------------------------------

#         # Perform 10-Fold Cross Validation

#         # Store the function's output as a variable
#         scores = af.get10FoldCVScore(X,y)

#         featBalancedAcc = scores.mean()
#         stdDev = scores.std()

#         df3 = df3.append({'feature': feature, 'iteration': iteration, 'featBalancedAcc': featBalancedAcc, 'stdDev': stdDev}, ignore_index=True)

#         print('Feature ', str(feature), ', Iteration ', str(iteration), '... Stored.')

# outFileName = 'featBalancedAcc_' + str(dataset) + '_shuff.txt'

# df3.to_csv(outFileName, index=False)
