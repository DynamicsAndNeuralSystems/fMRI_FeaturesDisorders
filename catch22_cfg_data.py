import glob
import scipy.io as sio
import numpy as np
import catch22

# Get all the feature names from the catch22 module in one list,
# excepting special and whole module attributes
features = dir(catch22)
features = [feat for feat in features if not '__' in feat]
del features[22:24] # last 2 attributes refer to the whole module.

#-------------------------------------------------------------------------------
def generatePathNames(path):
    ''' Returns sorted paths to all files with the extension '.mat' within the
    given folder.'''

    filePaths = sorted(glob.glob(path + '*.mat'))
    return filePaths
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
def fetchSubjectData(path,element):
    ''' Returns subject data (ts x roi ndarray) read from
    specified cell of 'roiTS' within the cfg data file '''

    matFile = sio.loadmat(path)
    subjectData = matFile['cfg'][0][0]['roiTS'][0][element-1]
    return subjectData
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
def catch22OneSeries(regionData):
    ''' Returns a catch-22 feature vector from a time-series '''
    global features

    featureVector = []
    for f in features:
        featureFunction = getattr(catch22,f)
        featureVector.append(featureFunction(regionData))

    return featureVector
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
def catch22OneRegion(filePaths,element,roi):
    ''' Returns a features x subjects matrix (ndarray) from the data of one
    region (roi) in all subjects. '''

    i = 0

    for path in filePaths:
        if i == 0:
            featureMatrix = np.zeros(shape=(len(filePaths),22))

        #fetch the subject's ts x roi array.
        subjectData = fetchSubjectData(path,element)
        # Fetch one column (i.e. one region of a subject)
        regionData = subjectData[:,roi].tolist()

        featureMatrix[i,:] = catch22OneSeries(regionData)
        i+=1

    return featureMatrix
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
def catch22AllRegions(path,element):
    '''Returns a features x subjects matrix (ndarray) with different regions
    simply one after the other in the matrix.
    I.e. to traverse from one region to another in the returned matrix,
    steps of size [no. of subjects] must be taken. '''

    # Generate a list of paths to each cfg data file
    filePaths = generatePathNames(path)

    # Get the shape of first subject's ts x roi array
    [timeRows, regionCols] = fetchSubjectData(filePaths[0],element).shape

    #Build feature x observation matrices, vertically stacking each region.
    for roi in range(regionCols):
        if roi == 0:
            allFeatMats = catch22OneRegion(filePaths,element,roi)

        else:
            featMat = catch22OneRegion(filePaths,element,roi)
            allFeatMats = np.vstack([allFeatMats,featMat])
    return allFeatMats
#-------------------------------------------------------------------------------

# Select a dataset:
dataset = 'COBRE'

# Assign the path to the dataset
path = '/Users/preethompal/Dropbox (Sydney Uni Student)/'+dataset+'/cfgData/'

#Select which brain parcellation will be analysed. See roiTS cells in cfgData files.
element = 1 #1 = first cell.

catch22CFGData = catch22AllRegions(path, element)

# Save file.
np.savetxt('c22' + '_roiTS' + str(element) + '_' + dataset + '.txt',catch22CFGData,delimiter=',')

'''NB File format: each line is a feature vector. Lines from the same
region across every subject are grouped together, i.e. First [no. of subjects]
lines belong to one region, then the next region etc. '''
