from selections import Selections
import after_catch_analysis_pack as acap
import pandas as pd

#------------------NB Make selections in selections.py before running-----------
s = Selections()
#Select framewise displacement threshold
threshold_fd = 0.5
# Select true to show distributions of t-statistics of every region for each feature. (I.e. Shows extent that sczs and ctrls have non-identical means)
showTValueHistograms = True
# Select true to graphically compare feature performance in both datasets. [Doesn't matter which dataset is chosen in selections.py]
showJointPlot = False
# Select true to show the estimated combined p values from both datasets using permutation testing.
showCombinedPValues = False # NB: Expensive, recommended to supress other functions if doing this.
usedSavedData = True #Used saved permutation testing results. Set to true this after running above once.
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
elif dataset == 'COBRE':
    labelColumn = acap.readLabelColumn(participants=participants,dataset='COBRE',indices=subjIndicesBelowThresh)

if showTValueHistograms:
    acap.tValueHistograms(dataset, featCount, featNames, roiCount, subjCount, c22Data, subjIndicesBelowThresh, labelColumn)

if showJointPlot or showCombinedPValues:
    if dataset == 'UCLA':
        alt = Selections(dataset='COBRE')
    elif dataset == 'COBRE':
        alt = Selections(dataset='UCLA')

    altFilePathsBelowThresh, altSubjIndicesBelowThresh = acap.removePathsAboveThresh(alt.fdAvgs, threshold_fd, alt.filePathsAll)

    if alt.dataset=='COBRE':
        altLabelColumn = acap.readLabelColumn(participants=alt.participants,dataset='COBRE',indices=altSubjIndicesBelowThresh)
    elif alt.dataset=='UCLA':
        altLabelColumn = acap.readLabelColumn(filePaths=altFilePathsBelowThresh,dataset='UCLA')

    featAccs = acap.accuracyOfFeatures(c22Data, roiCount, subjCount, featNames, subjIndicesBelowThresh, labelColumn, dispFigs, returnDF=True)
    featAccsAlt = acap.accuracyOfFeatures(alt.c22Data, alt.roiCount, alt.subjCount, alt.featNames, altSubjIndicesBelowThresh, altLabelColumn, alt.dispFigs, returnDF=True)

    if showJointPlot:
        acap.featurePerformanceJointPlot(featAccs, featAccsAlt, alt.dataset)
    if showCombinedPValues:
        if not usedSavedData:
            randomLearnData = acap.kiloLabelShufflesAndLearns(labelColumn, c22Data, subjIndicesBelowThresh, roiCount, subjCount, featNames)
            outFileName = "randomLearnData_"+dataset+".txt"
            randomLearnData.to_csv(outFileName, index=False)
            print("Your file "+outFileName+" has been saved in the current directory.")
            randomLearnData = 0
            altRandomLearnData = acap.kiloLabelShufflesAndLearns(altLabelColumn, alt.c22Data, altSubjIndicesBelowThresh, alt.roiCount, alt.subjCount, alt.featNames)
            outFileName = "randomLearnData_"+alt.dataset+".txt"
            altRandomLearnData.to_csv(outFileName, index=False)
            altRandomLearnData = 0
            print("Your file "+outFileName+" has been saved in the current directory.")
        randomLearnDataUCLA = pd.read_csv("randomLearnData_UCLA.txt")
        randomLearnDataCOBRE = pd.read_csv("randomLearnData_COBRE.txt")
        if dataset =='UCLA':
            featAccsUCLA = featAccs
            featAccsCOBRE = featAccsAlt
        else:
            featAccsUCLA = featAccsAlt
            featAccsCOBRE = featAccs
        # acap.featAccNullDistributionsPlot(randomLearnDataUCLA, randomLearnDataCOBRE, featAccsUCLA, featAccsCOBRE, featCount)
        acap.featureAccuracyPVals(featAccsUCLA,featAccsCOBRE,randomLearnDataUCLA,randomLearnDataCOBRE,featCount)
