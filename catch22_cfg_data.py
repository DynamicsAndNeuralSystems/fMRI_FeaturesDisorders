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
def retrieveExpData(matFile):
    subjectIDs = matFile['subject_list'][0]
    indices = []
    labels = []
    i = 0
    for subjectID in subjectIDs:
        if 'sub-1' in subjectID[0]:
            indices.append(i)
            labels.append(0)
        if 'sub-5' in subjectID[0]:
            indices.append(i)
            labels.append(1)
        i +=1
    allSubjectsData = matFile['time_series']
    sczsAndCtrls = allSubjectsData[:,:,indices,:]
    return sczsAndCtrls, labels


#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
def catch22FromSingleMatFile(path,element):
    element -= 1
    matFile = sio.loadmat(path)
    allSubjectsData, labels = retrieveExpData(matFile)
    numOfTimePoints, numOfRegions, numOfSubjects, numOfNoiseOptions = allSubjectsData.shape
    for roi in range(numOfRegions):
        regionSlice = allSubjectsData[:,roi,:,element]
        featureMatrix = np.zeros(shape=(numOfSubjects,22))
        for subject in range(numOfSubjects):
            subjectTimeSeries = regionSlice[:,subject].tolist()
            featureMatrix[subject,:]= catch22OneSeries(subjectTimeSeries)
        if roi == 0:
            allFeatMats = featureMatrix
        else:
            allFeatMats = np.vstack([allFeatMats,featureMatrix])
    return allFeatMats, labels
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
def fetchSubjectData(path,element):
    ''' Returns subject data (ts x roi ndarray) read from
    specified cell of 'roiTS' within a cfg data file '''

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
def catch22AllRegions(path,element,isSingleMatFile, dataOptionName):
    '''Returns a features x subjects matrix (ndarray) with different regions
    simply one after the other in the matrix.
    I.e. to traverse from one region to another in the returned matrix,
    steps of size [no. of subjects] must be taken. '''

    if not isSingleMatFile:
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

    if isSingleMatFile:
        allFeatMats, labels = catch22FromSingleMatFile(path,element)
        labelsFileName = 'labels_'+dataOptionName+str(element)+'_'+dataset+'.txt'
        np.savetxt(labelsFileName,labels,delimiter=',')
        print("Your file "+labelsFileName+" has been saved in the current directory.")

    return allFeatMats
#-------------------------------------------------------------------------------

# Select a dataset:
dataset = 'UCLA'
# Select True if all subjects' time series exist in a single mat file.
isSingleMatFile = True

if isSingleMatFile: #If using a single mat file select the path to the file here.
    path = '/Users/preethompal/Dropbox (Sydney Uni Student)/UCLA/UCLA_time_series_four_groups.mat'
else:
    # Assign the path to the dataset folder
    path = '/Users/preethompal/Dropbox (Sydney Uni Student)/'+dataset+'/cfgData/'

#Select which brain parcellation will be analysed. See 6 options of roiTS cells in cfgData files.
# If using single mat file (e.g. new DiCER data), use this to select which noise
# processing option is used. See shape of time_series (3 options)
element = 3 #1 = first option/cell.

dataOptionName = 'roiTS'
if isSingleMatFile:
    dataOptionName = 'procMeth'
catch22CFGData = catch22AllRegions(path, element, isSingleMatFile, dataOptionName)

# Save file.
saveFileName = 'c22_' + dataOptionName + str(element) + '_' + dataset + '.txt'
np.savetxt(saveFileName,catch22CFGData,delimiter=',')
print("Your file "+saveFileName+" has been saved in the current directory.")

'''NB File format: each line is a catch22 feature vector. Lines from the same
region across every subject are grouped together, i.e. First [no. of subjects]
lines belong to one region, then the next region etc. '''
