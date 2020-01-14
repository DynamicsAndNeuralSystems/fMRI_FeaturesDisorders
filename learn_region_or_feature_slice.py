import after_catch_analysis_pack as acap
from selections import Selections
'''Use this file to test classification accuracy region by region or feature by
feature; selecting each region or feature manually.'''
#------------------NB Make selections in selections.py before running-----------
s = Selections()

#Select the kind of slice you'd like to make; 'region' or 'feature'
sliceSelection = 'region'
# If looking at regions, choose which roi to analyse. Possible selections are 0 to [no. of regions-1]
roi = 0
# If looking at features, choose which feature to analyse. Possible selections are 0 to 21
feature = 0
# Select a threshold fd you'd like to use for this slice
threshold_fd = 0.5
# Set to true if you'd like to just do all region and feature slices at once (will show less detail).
doAll = True
#-------------------------------------------------------------------------------
dataset = s.dataset
featNames = s.featNames
participants = s.participants
filePathsAll = s.filePathsAll
subjCount = s.subjCount
c22Data = s.c22Data
roiCount = s.roiCount
featCount = s.featCount
fdAvgs = s.fdAvgs #Average framewise displacements for all subjects.
dispFigs = s.dispFigs

filePathsBelowThresh, subjIndicesBelowThresh = acap.removePathsAboveThresh(fdAvgs, threshold_fd, filePathsAll)
if dataset == 'UCLA':
    labelColumn = acap.readLabelColumn(filePaths=filePathsBelowThresh,dataset='UCLA')
if dataset == 'COBRE':
    labelColumn = acap.readLabelColumn(participants=participants,dataset='COBRE',indices=subjIndicesBelowThresh)

if doAll:
    acap.accuracyOfRegions(roiCount, c22Data, subjIndicesBelowThresh, labelColumn, dispFigs=True)
    acap.accuracyOfFeatures(c22Data, roiCount, subjCount, featNames, subjIndicesBelowThresh, labelColumn, dispFigs=True)
else:
    try:
        if sliceSelection == 'region':
            print("Doing region slice.")
            acap.regByRegAnalysis(roi, c22Data, labelColumn, roiCount, subjIndicesBelowThresh)
        elif sliceSelection == 'feature':
            print("Doing feature slice.")
            featureName = featNames[feature]
            acap.featByFeatAnalysis(featureName, c22Data, labelColumn, roiCount, subjCount, subjIndicesBelowThresh)
        else:
            print('Invalid slice selection.')
    except IndexError as e:
        print('Index error. Did you select a non-existing feature or region?')
