import pandas as pd
import numpy as np
from numpy import genfromtxt
import glob
import after_catch_analysis_pack as acap
from sklearn.model_selection import StratifiedKFold
from copy import copy
import matplotlib.pyplot as plt
import scipy.io as sio
# from catch22_cfg_data import catch22OneSeries

filePath = '/Users/preethompal/Dropbox (Sydney Uni Student)/COBRE/movementData/fdAvgs_COBRE.txt'
fName = 'c22_procMeth3_UCLA.txt'
c22Data = genfromtxt(fName, delimiter=',')
[subjRows, featCols] = c22Data.shape
featureListFileName = 'PythonFeatureList.txt' # This text file contains the names of the 22 features.
featNames = [lines.rstrip('\n') for lines in open(featureListFileName)]

def addMultiIndexTest(c22Data):
    subjCount = 171
    # Convert c22Data to a data frame with a MultiIndex. Return number of regions, subjects and features
    c22Data, roiCount, featCount = acap.addMultiIndex(c22Data,featNames, subjCount)
    return c22Data, roiCount, subjCount

def labelColumnTest():
    csvPath = '/Users/preethompal/Dropbox (Sydney Uni Student)/COBRE/participants.csv'
    indices2Keep = [0, 1, 2, 8, 10]
    participants = pd.read_csv(csvPath,header=0);
    participants.replace(to_replace=1,value=0,inplace=True)
    participants.replace(to_replace=2,value=1,inplace=True)
    targetCol = participants['Diagnosis'][indices2Keep]
    targetCol = np.asarray(targetCol,dtype=np.int).reshape(len(targetCol),1)
    print(targetCol[2])

def setIndexTest():
    df = pd.DataFrame({'month': [1, 4, 7, 10], 'year': [2012, 2014, 2013, 2014], 'sale': [55, 40, 84, 31]})
    df = df.set_index('month')
    print(df)

def featSliceTest(c22Data, featureName, roiCount, subjCount):
    subjIndicesBelowThresh = [0,1,2,3,4,5]
    roiIndices = list(range(0,roiCount))
    subjIndices = list(range(0,subjCount))
    featSlice = c22Data.loc[:,featureName].values.reshape(roiCount,subjCount).transpose()
    featSlice = pd.DataFrame(data=featSlice, index=subjIndices, columns=roiIndices).rename_axis('ROI', axis=1)
    featSlice.index.name = 'Subject'
    featSlice = featSlice.loc[subjIndicesBelowThresh,:]
    ctrlIndices = [0,1]
    print(featSlice.iloc[ctrlIndices,0])

def skfTest():
    skf = StratifiedKFold(n_splits=2)
    X = [[1,2,3],[3,4,6],[5,6,7], [6,4,9], [1,1,2], [2,8,7]]
    y = [0, 0, 1, 1, 0, 1]
    for train_index, test_index in (skf.split(X, y)):
        print(type(train_index), test_index)

def whereTest():
    a = np.array([[1],[0],[1],[0],[1]])
    print(a[0,0])

def copyTest():
    a = ["yo", "hi", "wassup"]
    b = copy(a)
    del b[2]
    print("a:", a)
    print("b:", b)

def plotAccuraciesVaryingWithFDTest():
    learnDataFilePath = 'learnData_roiTS1_'+'COBRE'+'_varyingFDThresh.txt'
    learnData = pd.read_csv(learnDataFilePath)
    avgROIError = learnData['AverageErrorUsingRegions']
    print(avgROIError)
    avgFeatError = learnData['AverageErrorUsingFeatures']
    fds = learnData['fd']
    sczToCtrlRatios = learnData['SCZ:Control']
    avgROIAccuracy = learnData['AverageAccuracyUsingRegions']
    avgFeatAccuracy = learnData['AverageAccuracyUsingFeatures']
    fig, ax1 = plt.subplots()
    ax1.plot(fds, avgFeatAccuracy,'b-')
    ax1.plot(fds, avgROIAccuracy, 'r-')
    plt.xlim(max(fds)+0.01, min(fds)-0.01)
    ax1.set_xlabel('fd (cm)')
    ax1.set_ylabel('Balanced Acc (%)')

    ax2 = ax1.twinx()
    ax2.plot(fds, sczToCtrlRatios,'g-')
    ax2.set_ylabel('SCZ:Control', color='g')
    ax2.tick_params('y', colors='g')
    plt.grid(b=None)

    ax1.legend(['Avg Feat Acc','Avg ROI Acc'],loc=3)
    ax2.legend(['SCZ:Control'],loc=2)
    fig.tight_layout()
    plt.show()

# def readingMatFileTest():
#     path = '/Users/preethompal/Dropbox (Sydney Uni Student)/UCLA/UCLA_time_series_four_groups.mat'
#     element = 3
#     element -= 1
#     matFile = sio.loadmat(path)
#     allSubjectsData = matFile['time_series']
#     numOfTimePoints, numOfRegions, numOfSubjects, numOfNoiseOptions = allSubjectsData.shape
#     for roi in range(numOfRegions):
#         regionSlice = allSubjectsData[:,roi,:,element]
#         featureMatrix = np.zeros(shape=(numOfSubjects,22))
#         for subject in range(numOfSubjects):
#             subjectTimeSeries = regionSlice[:,subject].tolist()
#             featureMatrix[subject,:]= catch22OneSeries(subjectTimeSeries)
#         if roi == 0:
#             allFeatMats = featureMatrix
#         else:
#             allFeatMats = np.vstack([allFeatMats,featureMatrix])

def UCLAParticipantsTest():
    csvPath = '/Users/preethompal/Dropbox (Sydney Uni Student)/UCLA/participants.csv'
    participants = pd.read_csv(csvPath,header=0);
    print(len(participants))
def UCLAFilePathsTest():
    folderPath = '/Users/preethompal/Dropbox (Sydney Uni Student)/UCLA/cfgData/'
    filePaths = sorted(glob.glob(folderPath + '*.mat'))
    print(len(filePaths))
def subjSliceTest(c22Data, subjCount, roiCount):
    for i in range(subjCount):
        subjSlice = c22Data.iloc[c22Data.index.get_level_values(1)==i]
        for j in range(roiCount):
            regionInSubject = np.asarray(subjSlice.iloc[j])
            print(regionInSubject)
            break
        break

c22Data, roiCount, subjCount = addMultiIndexTest(c22Data)
subjSliceTest(c22Data, subjCount, roiCount)

# featSliceTest(c22Data, featNames[0], roiCount, subjCount)
