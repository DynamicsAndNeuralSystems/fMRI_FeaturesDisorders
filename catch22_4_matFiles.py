import glob
import scipy.io as sio
import catch22
import numpy as np

path = '/Users/AV/Dropbox/UCLA/'

# This will arrange the assigned file type in ALPHABETICAL ORDER
TS_path_names = sorted(glob.glob(path + '*.mat'))

#-------------------------------------------------------------------------------
''' This function extracts data from any element of 'roiTS' within the cfg
    structure in the .mat file '''

def gettsData(path,element):
    matFile = sio.loadmat(path)
    tsData = matFile['cfg'][0][0]['roiTS'][0][element]
    return tsData
#-------------------------------------------------------------------------------
''' This function returns a catch-22 feature vector from an input time-series
    list '''

def giveMeFeatureVector(tsData):

    features = dir(catch22)
    features = [item for item in features if not '__' in item]

    featureVector = []
    for testFun in features:
        featureFun = getattr(catch22,testFun)
        featureVector.append(featureFun(tsData))

    return featureVector
#-------------------------------------------------------------------------------
''' This function takes any element of the 'roiTS' cell from each of the
    subject data files and an input ROI (a specified column of the matrix),
    converts these columns into a list, processes them using the catch22 module
    and stores them in a variable that is returned, smallFeatMat '''

def giveMeSmallFeatMat(element,roi):

    # Initialise indices
    i = 0

    for path in TS_path_names:

        # Retrieve and store data in a temp variable
        tsData = gettsData(path,element)

        # Need to initialise the smallFeatMat for the first iteration of the loop
        if i == 0:

            smallFeatMat = np.zeros(shape=(len(TS_path_names),22))

        # Get a specific column from the ts x roi data matrix for all the files
        tsDataCol = tsData[:,roi]

        # To be recognised by the catch22 modules, the input array needs to be a list
        tsDataCol = tsDataCol.tolist()
        smallFeatMat[i,:] = giveMeFeatureVector(tsDataCol)

        # Increment index
        i += 1

    return smallFeatMat
#-------------------------------------------------------------------------------
''' This function takes the index of an element in the cell 'roiTS' as an input
    and returns a matrix composed of rows of 1 x 22 feature vectors for every
    subject for each ROI in the specified element
    Note: The number of ROI varies from element to element '''

def giveMeBigFeatMat(element):

    # Assuming each of the data files have an indentical structure, use the first
    # element's data as a template to initialise a few variables
    [rows, cols] = gettsData(TS_path_names[0],element).shape

    for n in range(cols):

        if n == 0:
            bigFeatMat = giveMeSmallFeatMat(element,n)

        else:
            smallFeatMat = giveMeSmallFeatMat(element,n)
            bigFeatMat = np.vstack([bigFeatMat,smallFeatMat])

    return bigFeatMat
#-------------------------------------------------------------------------------

# Testing...
x = giveMeBigFeatMat(2)

# np.savetxt('element2.txt',x,delimiter=',')
