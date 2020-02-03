import after_catch_analysis_pack as acap
import numpy as np
import pandas as pd
from selections import Selections
'''Use this file to test classification accuracy region by region or feature by
feature; selecting each region or feature manually.'''
#------------------NB Make selections in selections.py before running-----------
s = Selections()

#Select the kind of slice you'd like to make; 'region' or 'feature'
sliceSelection = 'feature'
# If looking at regions, choose which roi to analyse. Possible selections are 0 to [no. of regions-1]
roi = 0
# If looking at features, choose which feature to analyse. Possible selections are 0 to 21
feature = 16
# Select a threshold fd you'd like to use for this slice
threshold_fd = 0.5
# Set to true if you'd like to just do all region and feature slices at once (will show all accuracies in histogram).
doAll = True
doFeaturesAndRegionsAtOnce = False
useSavedAccuracies = True
#-------------------------------------------------------------------------------
dataset = s.dataset
dataOptionName = s.dataOptionName
featNames = s.featNames
participants = s.participants
filePathsAll = s.filePathsAll
subjCount = s.subjCount
problemSubjectInds = s.problemSubjectInds
c22Data = s.c22Data
roiCount = s.roiCount
featCount = s.featCount
fdAvgs = s.fdAvgs #Average framewise displacements for all subjects.
dispFigs = s.dispFigs
# filePathsBelowThresh, subjIndicesBelowThresh = acap.removePathsAboveThresh(fdAvgs, threshold_fd, filePathsAll)

if 'procMeth' in dataOptionName:
    print("Reading from "+dataOptionName)
    labelColumn = np.genfromtxt('labels_'+dataOptionName+'_UCLA.txt').astype(int)
    labelColumn = np.reshape(labelColumn,(len(labelColumn), 1))
    c22Data, labelColumn, subjCount = acap.killNaNs(c22Data, labelColumn, problemSubjectInds)
    subjIndicesBelowThresh = [i for i in range(subjCount)]
    print("Subj Count:", subjCount, "\n")

    dicer2 = Selections("UCLA", "procMeth2")
    print("Reading from "+dicer2.dataOptionName)
    labelColumnDicer2 = np.genfromtxt('labels_'+dicer2.dataOptionName+'_UCLA.txt').astype(int)
    labelColumnDicer2 = np.reshape(labelColumnDicer2,(len(labelColumnDicer2), 1))
    dicer2.c22Data, labelColumnDicer2, dicer2.subjCount = acap.killNaNs(dicer2.c22Data, labelColumnDicer2, dicer2.problemSubjectInds)
    subjIndicesBelowThreshDicer2 = [i for i in range(dicer2.subjCount)]
    print("Subj Count:", dicer2.subjCount, "\n")

    dicer3 = Selections("UCLA", "procMeth3")
    print("Reading from "+dicer3.dataOptionName)
    labelColumnDicer3 = np.genfromtxt('labels_'+dicer3.dataOptionName+'_UCLA.txt').astype(int)
    labelColumnDicer3 = np.reshape(labelColumnDicer3,(len(labelColumnDicer3), 1))
    dicer3.c22Data, labelColumnDicer3, dicer3.subjCount = acap.killNaNs(dicer3.c22Data, labelColumnDicer3, dicer3.problemSubjectInds)
    subjIndicesBelowThreshDicer3 = [i for i in range(dicer3.subjCount)]
    print("Subj Count:", dicer3.subjCount, "\n")


elif dataset == 'UCLA':
    labelColumn = acap.readLabelColumn(filePaths=filePathsBelowThresh,dataset='UCLA')
if dataset == 'COBRE':
    labelColumn = acap.readLabelColumn(participants=participants,dataset='COBRE',indices=subjIndicesBelowThresh)

if doAll:
    if not useSavedAccuracies:
        print("\nLearning from "+dataOptionName+'\n')
        roiAccs = acap.accuracyOfRegions(roiCount, c22Data, subjIndicesBelowThresh, labelColumn, dispFigs=False, returnDF=True)
        outFileName = 'roiAccuracies_procMeth1_UCLA.txt'
        roiAccs.to_csv(outFileName, index=False, header=True)
        roiAccs = 0

        featAccs = acap.accuracyOfFeatures(c22Data, roiCount, subjCount, featNames, subjIndicesBelowThresh, labelColumn, dispFigs=False, returnDF=True)
        outFileName = 'featAccuracies_procMeth1_UCLA.txt'
        featAccs.to_csv(outFileName, index=False, header=True)
        featAccs = 0

        print("\nLearning from "+dicer2.dataOptionName+'\n')
        roiAccsDicer2 = acap.accuracyOfRegions(dicer2.roiCount, dicer2.c22Data, subjIndicesBelowThreshDicer2, labelColumnDicer2, dispFigs=False, returnDF=True)
        outFileName = 'roiAccuracies_procMeth2_UCLA.txt'
        roiAccsDicer2.to_csv(outFileName, index=False, header=True)
        roiAccsDicer2 = 0

        featAccsDicer2 = acap.accuracyOfFeatures(dicer2.c22Data, dicer2.roiCount, dicer2.subjCount, featNames,
                                                subjIndicesBelowThreshDicer2, labelColumnDicer2, dispFigs=False, returnDF=True)
        outFileName = 'featAccuracies_procMeth2_UCLA.txt'
        featAccsDicer2.to_csv(outFileName, index=False, header=True)
        featAccsDicer2 = 0

        print("\nLearning from "+dicer3.dataOptionName+'\n')
        roiAccsDicer3 = acap.accuracyOfRegions(dicer3.roiCount, dicer3.c22Data, subjIndicesBelowThreshDicer3, labelColumnDicer3, dispFigs=False, returnDF=True)
        outFileName = 'roiAccuracies_procMeth3_UCLA.txt'
        roiAccsDicer3.to_csv(outFileName, index=False, header=True)
        roiAccsDicer3 = 0

        featAccsDicer3 = acap.accuracyOfFeatures(dicer3.c22Data, dicer3.roiCount, dicer3.subjCount, featNames,
                                                subjIndicesBelowThreshDicer3, labelColumnDicer3, dispFigs=False, returnDF=True)
        outFileName = 'featAccuracies_procMeth3_UCLA.txt'
        featAccsDicer3.to_csv(outFileName, index=False, header=True)
        featAccsDicer3 = 0

    featAccs = pd.read_csv('featAccuracies_procMeth1_UCLA.txt')
    featAccsDicer2 = pd.read_csv('featAccuracies_procMeth2_UCLA.txt')
    featAccsDicer3 = pd.read_csv('featAccuracies_procMeth3_UCLA.txt')

    roiAccs = pd.read_csv('roiAccuracies_procMeth1_UCLA.txt')
    roiAccsDicer2 = pd.read_csv('roiAccuracies_procMeth2_UCLA.txt')
    roiAccsDicer3 = pd.read_csv('roiAccuracies_procMeth3_UCLA.txt')

    acap.distributionsOverlayPlot([roiAccs, roiAccsDicer2, roiAccsDicer3], '% Accuracy', 'Regions')
    acap.distributionsOverlayPlot([featAccs, featAccsDicer2, featAccsDicer3], '% Accuracy', 'Features')

elif doFeaturesAndRegionsAtOnce:
    acap.accuracyOfRegionsAndFeatures(c22Data, labelColumn, roiCount, subjCount, featCount)
    acap.accuracyOfRegionsAndFeatures(dicer2.c22Data, labelColumnDicer2, dicer2.roiCount, dicer2.subjCount, dicer2.featCount)
    acap.accuracyOfRegionsAndFeatures(dicer3.c22Data, labelColumnDicer3, dicer3.roiCount, dicer3.subjCount, dicer3.featCount)
else:
    try:
        if sliceSelection == 'region':
            print("Doing region slice.")
            acap.regByRegAnalysis(roi, c22Data, labelColumn, roiCount, subjIndicesBelowThresh)
        elif sliceSelection == 'feature':
            print("Doing feature slice.")
            featureName = featNames[feature]
            print("\nLearning from "+dataOptionName+'\n')
            acap.featByFeatAnalysis(featureName, c22Data, labelColumn, roiCount, subjCount, subjIndicesBelowThresh)
            print("\nLearning from "+dicer2.dataOptionName+'\n')
            acap.featByFeatAnalysis(featureName, dicer2.c22Data, labelColumnDicer2, dicer2.roiCount, dicer2.subjCount, subjIndicesBelowThreshDicer2)
            print("\nLearning from "+dicer3.dataOptionName+'\n')
            acap.featByFeatAnalysis(featureName, dicer3.c22Data, labelColumnDicer3, dicer3.roiCount, dicer3.subjCount, subjIndicesBelowThreshDicer3)
        else:
            print('Invalid slice selection.')
    except IndexError as e:
        print('Index error. Did you select a non-existing feature or region?')
