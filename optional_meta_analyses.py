from selections import Selections
import after_catch_analysis_pack as acap
import pandas as pd
import numpy as np
#------------------NB Make selections in selections.py before running-----------
s = Selections()
#Select framewise displacement threshold
threshold_fd = 0.5
# Select true to show distribution plots of where NaNs exist
showNaNPlots = False
# Select true to show distributions of t-statistics of every region for each feature. (I.e. Shows extent that sczs and ctrls have non-identical means)
showTValueHistograms = False
# Select true to graphically compare feature performance in both datasets. [Doesn't matter which dataset is chosen in selections.py]
showJointPlot = False
# Select true to show the estimated p values from both datasets' accuracies using permutation testing.
showIndividualFeatureSignificance = False # NB: Expensive if getting random permutations first time.
showIndividualRegionSignificance = False
showJointFeatureSignificance = False
showJointRegionSignificance = False
showJointFeaturesAndRegionsSignificance = True

useSavedRandomLearnData = False #Used saved permutation testing results. Set to true this after running above once.
useSavedAccuracies = True
#-------------------------------------------------------------------------------

dataset = s.dataset
dataOptionName = s.dataOptionName
featNames = s.featNames
participants = s.participants
filePathsAll = s.filePathsAll
subjCount = s.subjCount
c22Data = s.c22Data
roiCount = s.roiCount
featCount = s.featCount
fdAvgs = s.fdAvgs #Average framewise displacements for all subjects.
dispFigs = s.dispFigs
problemSubjectInds = s.problemSubjectInds
# filePathsBelowThresh, subjIndicesBelowThresh = acap.removePathsAboveThresh(fdAvgs, threshold_fd, filePathsAll)

if 'procMeth' in dataOptionName:
    print("Reading from "+dataOptionName)
    labelColumn = np.genfromtxt('labels_'+dataOptionName+'_UCLA.txt').astype(int)
    labelColumn = np.reshape(labelColumn,(len(labelColumn), 1))
    c22Data, labelColumn, subjCount = acap.killNaNs(c22Data, labelColumn, problemSubjectInds)
    subjIndicesBelowThresh = [i for i in range(subjCount)]


    dicer2 = Selections("UCLA", "procMeth2")
    print("Reading from "+dicer2.dataOptionName)
    labelColumnDicer2 = np.genfromtxt('labels_'+dicer2.dataOptionName+'_UCLA.txt').astype(int)
    labelColumnDicer2 = np.reshape(labelColumnDicer2,(len(labelColumnDicer2), 1))
    dicer2.c22Data, labelColumnDicer2, dicer2.subjCount = acap.killNaNs(dicer2.c22Data, labelColumnDicer2, dicer2.problemSubjectInds)
    subjIndicesBelowThreshDicer2 = [i for i in range(dicer2.subjCount)]


    dicer3 = Selections("UCLA", "procMeth3")
    print("Reading from "+dicer3.dataOptionName)
    labelColumnDicer3 = np.genfromtxt('labels_'+dicer3.dataOptionName+'_UCLA.txt').astype(int)
    labelColumnDicer3 = np.reshape(labelColumnDicer3,(len(labelColumnDicer3), 1))
    dicer3.c22Data, labelColumnDicer3, dicer3.subjCount = acap.killNaNs(dicer3.c22Data, labelColumnDicer3, dicer3.problemSubjectInds)
    subjIndicesBelowThreshDicer3 = [i for i in range(dicer3.subjCount)]

    # roiInds = [37,45,78]
    # weirdVarSlice = np.asarray(c22Data.loc[roiInds,:]['CO_f1ecac'])

if showNaNPlots:
    acap.nanDistribution(c22Data, subjCount, roiCount, featCount,disp=True)

if showTValueHistograms:
    tValDFs = acap.tValueHistograms(dataset, featCount, featNames, roiCount, subjCount, c22Data, subjIndicesBelowThresh, labelColumn, returnDFs=True)
    tValDFsDicer2 = acap.tValueHistograms(dicer2.dataset, dicer2.featCount, featNames, dicer2.roiCount, dicer2.subjCount, dicer2.c22Data, subjIndicesBelowThreshDicer2, labelColumnDicer2, returnDFs=True)
    tValDFsDicer3 = acap.tValueHistograms(dicer3.dataset, dicer3.featCount, featNames, dicer3.roiCount, dicer3.subjCount, dicer3.c22Data, subjIndicesBelowThreshDicer3, labelColumnDicer3, returnDFs=True)

    acap.tValueHistogramsOverlay([tValDFs,tValDFsDicer2,tValDFsDicer3], featCount)

if showJointPlot:
    acap.featurePerformanceJointPlot(featAccs, featAccsAlt, alt.dataset)

if showIndividualFeatureSignificance or showJointFeatureSignificance:
    if not useSavedAccuracies:
        featAccs = acap.accuracyOfFeatures(c22Data, roiCount, subjCount, featNames, subjIndicesBelowThresh, labelColumn, dispFigs, returnDF=True)
        outFileName = 'featAccuracies_procMeth1_UCLA.txt'
        featAccs.to_csv(outFileName, index=False, header=True)
        featAccs = 0

        featAccsDicer2 = acap.accuracyOfFeatures(dicer2.c22Data, dicer2.roiCount, dicer2.subjCount, featNames,
                                                subjIndicesBelowThreshDicer2, labelColumnDicer2, dispFigs=False, returnDF=True)
        outFileName = 'featAccuracies_procMeth2_UCLA.txt'
        featAccsDicer2.to_csv(outFileName, index=False, header=True)
        featAccsDicer2 = 0

        featAccsDicer3 = acap.accuracyOfFeatures(dicer3.c22Data, dicer3.roiCount, dicer3.subjCount, featNames,
                                                subjIndicesBelowThreshDicer3, labelColumnDicer3, dispFigs=False, returnDF=True)
        outFileName = 'featAccuracies_procMeth3_UCLA.txt'
        featAccsDicer3.to_csv(outFileName, index=False, header=True)
        featAccsDicer3 = 0


    if not useSavedRandomLearnData and showIndividualFeatureSignificance:
        randomLearnData = acap.kiloLabelShufflesAndLearnsFeatures(labelColumnDicer3, dicer3.c22Data, subjIndicesBelowThreshDicer3, dicer3.roiCount, dicer3.subjCount, featNames)
        outFileName = "randomLearnData_individualFeatures_procMeth3_UCLA.txt"
        randomLearnData.to_csv(outFileName, index=False, header=True)
        print("Your file "+outFileName+" has been saved in the current directory.")
        randomLearnData = 0

    if not useSavedRandomLearnData and showJointFeatureSignificance:
        randomLearnData = acap.kiloLabelShufflesAndLearnsFeaturesJoint(labelColumnDicer3, dicer3.c22Data, subjIndicesBelowThreshDicer3, dicer3.roiCount, dicer3.subjCount, featNames)
        outFileName = "randomLearnData_jointFeatures_procMeth3_UCLA.txt"
        randomLearnData.to_csv(outFileName, index=False, header=True)
        print("Your file "+outFileName+" has been saved in the current directory.")
        randomLearnData = 0

    featAccs = pd.read_csv('featAccuracies_procMeth1_UCLA.txt')
    featAccsDicer2 = pd.read_csv('featAccuracies_procMeth2_UCLA.txt')
    featAccsDicer3 = pd.read_csv('featAccuracies_procMeth3_UCLA.txt')

    if showIndividualFeatureSignificance:
        randomLearnDataDicerUCLA = pd.read_csv("randomLearnData_individualFeatures_procMeth3_UCLA.txt")
        acap.featAccNullDistributionsPlotTriple(randomLearnDataDicerUCLA, [featAccs, featAccsDicer2, featAccsDicer3], featCount)
        acap.featureAccuracyPValsTriple([featAccs, featAccsDicer2, featAccsDicer3], randomLearnDataDicerUCLA, featCount)

    if showJointFeatureSignificance:
        randomLearnDataDicerUCLA = pd.read_csv("randomLearnData_jointFeatures_procMeth3_UCLA.txt")
        meanFeatAcc = np.mean(np.asarray(featAccs['% Accuracy']))
        meanFeatAccDicer2 = np.mean(np.asarray(featAccsDicer2['% Accuracy']))
        meanFeatAccDicer3 = np.mean(np.asarray(featAccsDicer3['% Accuracy']))
        acap.jointAccNullDistributionPlot([meanFeatAcc, meanFeatAccDicer2, meanFeatAccDicer3], randomLearnDataDicerUCLA, 'Features')
        print(meanFeatAcc, meanFeatAccDicer2, meanFeatAccDicer3)
        acap.jointAccuracyPValTriple([meanFeatAcc, meanFeatAccDicer2, meanFeatAccDicer3], randomLearnDataDicerUCLA)

if showIndividualRegionSignificance or showJointRegionSignificance:
    if not useSavedAccuracies:
        roiAccs = acap.accuracyOfRegions(roiCount, c22Data, subjIndicesBelowThresh, labelColumn, dispFigs=False, returnDF=True)
        outFileName = 'roiAccuracies_procMeth1_UCLA.txt'
        roiAccs.to_csv(outFileName, index=False, header=True)
        roiAccs = 0

        roiAccsDicer2 = acap.accuracyOfRegions(dicer2.roiCount, dicer2.c22Data, subjIndicesBelowThreshDicer2, labelColumnDicer2, dispFigs=False, returnDF=True)
        outFileName = 'roiAccuracies_procMeth2_UCLA.txt'
        roiAccsDicer2.to_csv(outFileName, index=False, header=True)
        roiAccsDicer2 = 0

        roiAccsDicer3 = acap.accuracyOfRegions(dicer3.roiCount, dicer3.c22Data, subjIndicesBelowThreshDicer3, labelColumnDicer3, dispFigs=False, returnDF=True)
        outFileName = 'roiAccuracies_procMeth3_UCLA.txt'
        roiAccsDicer3.to_csv(outFileName, index=False, header=True)
        roiAccsDicer3 = 0

    if not useSavedRandomLearnData and showIndividualRegionSignificance:
        randomLearnData = acap.kiloLabelShufflesAndLearnsRegions(labelColumnDicer3, dicer3.c22Data, dicer3.roiCount)
        outFileName = "randomLearnData_individualRegions_procMeth3_UCLA.txt"
        randomLearnData.to_csv(outFileName, mode='a', index=False, header=False)
        print("Your file "+outFileName+" has been saved in the current directory.")
        randomLearnData = 0

    if not useSavedRandomLearnData and showJointRegionSignificance:
        randomLearnData = acap.kiloLabelShufflesAndLearnsRegionsJoint(labelColumnDicer3, dicer3.c22Data, dicer3.roiCount)
        outFileName = 'randomLearnData_jointRegions_procMeth3_UCLA.txt'
        randomLearnData.to_csv(outFileName, index=False, header=True)
        print("Your file "+outFileName+" has been saved in the current directory.")
        randomLearnData = 0

    roiAccs = pd.read_csv('roiAccuracies_procMeth1_UCLA.txt')
    roiAccsDicer2 = pd.read_csv('roiAccuracies_procMeth2_UCLA.txt')
    roiAccsDicer3 = pd.read_csv('roiAccuracies_procMeth3_UCLA.txt')

    if showIndividualRegionSignificance:
        randomLearnDataDicerUCLA = pd.read_csv("randomLearnData_individualRegions_procMeth3_UCLA.txt")
        acap.roiAccNullDistributionsPlotTriple(randomLearnDataDicerUCLA, [roiAccs, roiAccsDicer2, roiAccsDicer3], roiCount)
        acap.roiAccuracyPValsTriple(randomLearnDataDicerUCLA, [roiAccs, roiAccsDicer2, roiAccsDicer3], roiCount)

    if showJointRegionSignificance:
        randomLearnDataDicerUCLA = pd.read_csv("randomLearnData_jointRegions_procMeth3_UCLA.txt")
        meanRoiAcc = np.mean(np.asarray(roiAccs['% Accuracy']))
        meanRoiAccDicer2 = np.mean(np.asarray(roiAccsDicer2['% Accuracy']))
        meanRoiAccDicer3 = np.mean(np.asarray(roiAccsDicer3['% Accuracy']))
        # acap.jointAccNullDistributionPlot([meanRoiAcc, meanRoiAccDicer2, meanRoiAccDicer3], randomLearnDataDicerUCLA, 'Regions')
        acap.jointAccuracyPValTriple([meanRoiAcc, meanRoiAccDicer2, meanRoiAccDicer3], randomLearnDataDicerUCLA)

if showJointFeaturesAndRegionsSignificance:
    jointRegionsFeaturesAcc = acap.accuracyOfRegionsAndFeatures(c22Data, labelColumn, roiCount, subjCount, featCount, returnAcc=True)
    jointRegionsFeaturesAccDicer2 = acap.accuracyOfRegionsAndFeatures(dicer2.c22Data, labelColumnDicer2, dicer2.roiCount, dicer2.subjCount, dicer2.featCount, returnAcc=True)
    jointRegionsFeaturesAccDicer3 = acap.accuracyOfRegionsAndFeatures(dicer3.c22Data, labelColumnDicer3, dicer3.roiCount, dicer3.subjCount, dicer3.featCount, returnAcc=True)

    if not useSavedRandomLearnData:
        randomLearnData = acap.kiloLabelShufflesAndLearnsRegionsAndFeatures(dicer3.c22Data, labelColumnDicer3, dicer3.roiCount, dicer3.subjCount, dicer3.featCount)
        outFileName = 'randomLearnData_FeaturesAndRegions_procMeth3_UCLA.txt'
        randomLearnData.to_csv(outFileName, index=False,header=True)
        print("Your file "+outFileName+" has been saved in the current directory.")
        randomLearnData = 0

    randomLearnDataDicerUCLA = pd.read_csv('randomLearnData_FeaturesAndRegions_procMeth3_UCLA.txt')

    print(jointRegionsFeaturesAcc,jointRegionsFeaturesAccDicer2,jointRegionsFeaturesAccDicer3)
    # acap.jointAccNullDistributionPlot([jointRegionsFeaturesAcc,jointRegionsFeaturesAccDicer2,jointRegionsFeaturesAccDicer3], randomLearnDataDicerUCLA, 'Regions and Features')
    acap.jointAccuracyPValTriple([jointRegionsFeaturesAcc,jointRegionsFeaturesAccDicer2,jointRegionsFeaturesAccDicer3], randomLearnDataDicerUCLA)
