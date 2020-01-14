import pandas as pd
import numpy as np
import after_catch_analysis_pack as acap
from selections import Selections

'''Use this file to test average classification accuracy of both all regions and
all features, for each given threshold framewise displacement. This is
stored in a data frame and saved to csv. '''

'''Turn on dispFigs in selections.py to show a plot of accuracies of all regions,
then of all features for each threshold fd as they are calculated. Will slow down.'''

'''Once this program has been run once, you can view plots from the file saved
again by setting useSavedData to True'''
#------------------NB Make selections in selections.py before running-----------
s = Selections()
useSavedData = True
#-------------------------------------------------------------------------------

fdAvgs = s.fdAvgs #Average framewise displacements for all subjects.
dataset = s.dataset
if not useSavedData:
    dispFigs = s.dispFigs
    featNames = s.featNames
    participants = s.participants
    filePathsAll = s.filePathsAll
    subjCount = s.subjCount
    c22Data = s.c22Data
    roiCount = s.roiCount
    featCount = s.featCount
    threshold_fds = np.linspace(0.72,0.21,18) # Generate possible thresholds.

    df = pd.DataFrame({'fd': [], 'SCZ': [], 'Control': [], 'Total': [],
                        'SCZ:Control': [], 'AverageAccuracyUsingRegions': [], 'AverageErrorUsingRegions': [], 'AverageAccuracyUsingFeatures': [], 'AverageErrorUsingFeatures': []})
    for threshold_fd in threshold_fds:
        # Get file paths and indices to subjects that have an fd < threshold_fd
        filePathsBelowThresh, subjIndicesBelowThresh = acap.removePathsAboveThresh(fdAvgs, threshold_fd, filePathsAll)

        # Create the label column (ndarray) for those below this threshold
        if dataset == 'UCLA':
            labelColumn = acap.readLabelColumn(filePaths=filePathsBelowThresh,dataset='UCLA')
        if dataset == 'COBRE':
            labelColumn = acap.readLabelColumn(participants=participants,dataset='COBRE',indices=subjIndicesBelowThresh)
        ctrlCount, sczCount, totalCount, sczToCtrlRatio = acap.groupCounts(labelColumn)

        # Get average accuracy and error for a region to classify ctrl vs scz.
        # Can turn on figure display to show best regions for classification.
        avgROIAcc, avgROIError = acap.accuracyOfRegions(roiCount, c22Data, subjIndicesBelowThresh, labelColumn, dispFigs)
        # Get average accuracy and error for a c22 feature to classify ctrl vs scz.
        # Can turn on figure display to show best features for classification.
        avgFeatAcc, avgFeatError = acap.accuracyOfFeatures(c22Data, roiCount, subjCount, featNames, subjIndicesBelowThresh, labelColumn, dispFigs)

        df = df.append({'fd': threshold_fd, 'SCZ': sczCount, 'Control': ctrlCount,
                        'Total': totalCount, 'SCZ:Control': sczToCtrlRatio,
                        'AverageAccuracyUsingRegions': avgROIAcc,
                        'AverageErrorUsingRegions': avgROIError,
                        'AverageAccuracyUsingFeatures': avgFeatAcc,
                        'AverageErrorUsingFeatures': avgFeatError}, ignore_index=True)

    saveFileName = 'learnData_roiTS1_'+dataset+'_varyingFDThresh.txt'
    df.to_csv(saveFileName, index=False)
    print("Your file "+saveFileName+" has been saved in the current directory.")
    acap.accuraciesVaryingWithFDPlot(df)
    acap.showMeFDAcrossSubs(df,fdAvgs)

else:
    learnDataFilePath = 'learnData_roiTS1_'+dataset+'_varyingFDThresh.txt'
    learnData = pd.read_csv(learnDataFilePath)
    acap.accuraciesVaryingWithFDPlot(learnData)
    acap.showMeFDAcrossSubs(learnData,fdAvgs)
