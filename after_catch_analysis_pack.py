from copy import copy
import glob
import scipy.io as sio
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('ggplot')
import seaborn as sns
import math
from sklearn.svm import SVC
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.metrics import balanced_accuracy_score
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from scipy.stats import zscore
from scipy import stats
# from sm.stats import multitest
import statsmodels as sm
#-------------------------------------------------------------------------------
def removePathsAboveThresh(fdAvgs, threshold_fd, filePaths):
    ''' This function removes the entries within filePaths if those subjects had
    a framewise displacement higher than the threshold. '''

    indices2Del = np.where(fdAvgs > threshold_fd)[0]
    indices2Keep = np.where(fdAvgs <= threshold_fd)[0].tolist()
    tempFilePaths = copy(filePaths)
    for i in sorted(indices2Del, reverse=True):
        del tempFilePaths[i]
    return tempFilePaths, indices2Keep
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def nanDistribution(c22Data, subjCount, roiCount, featCount, disp=True):
    problemSubjectInds = []
    problemSubjectFeatureLocs = []
    problemSubjectAvgNumOfNans = []

    for i in range(subjCount):
        subjSlice = c22Data.iloc[c22Data.index.get_level_values('Subject')==i]
        numOfNansInEachRegion = np.zeros(roiCount)

        for j in range(roiCount):
            roiInSubject = np.asarray(subjSlice.iloc[j])
            featureLocOfNans = np.isnan(roiInSubject)
            numOfNansInEachRegion[j] = featureLocOfNans.sum()
            if j == 0:
                featureLocOfNansInEachRegion = featureLocOfNans
            else:
                featureLocOfNansInEachRegion = np.vstack((featureLocOfNansInEachRegion, featureLocOfNans))

        avgNumOfNansInSubject = numOfNansInEachRegion.mean()
        if avgNumOfNansInSubject > 0:
            featureLocOfNans = np.zeros(featCount)
            for j in range(featCount):
                featureLocOfNans[j] = featureLocOfNansInEachRegion[:,j].mean()
            problemSubjectInds.append(i)
            problemSubjectFeatureLocs.append(featureLocOfNans)
            problemSubjectAvgNumOfNans.append(avgNumOfNansInSubject)

    if disp:
        print("Subjects with at least one NaN somewhere:", problemSubjectInds , '\n')
        print("Showing Average Presence of NaNs in each feature across all regions (%) for these problem subjects...\n")
        df = pd.DataFrame(data=np.transpose(problemSubjectFeatureLocs), columns=[i for i in problemSubjectInds])
        df.rename_axis('Feature', axis='index',inplace=True)
        df.rename_axis('Subjects', axis='columns',inplace=True)
        df = df.apply(lambda x: x*100)
        print(df)
    return problemSubjectInds
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def killNaNs(c22Data, labelColumn, problemSubjectInds):
    print("Removing these problem subjects. Use showNanDistribution to see presence of NaNs.", problemSubjectInds)
    labelColumn = np.delete(labelColumn, problemSubjectInds)
    labelColumn = labelColumn.reshape((len(labelColumn),1))
    newSubjCount = len(labelColumn)

    c22Data = c22Data.drop(index=problemSubjectInds,level='Subject')
    newMultiIndex = []
    newMultiIndex.append(list(c22Data.index.levels[0])) #ROI indices
    newMultiIndex.append([i for i in range(newSubjCount)]) #Subj Indices
    c22Data.index = pd.MultiIndex.from_product(newMultiIndex, names=['ROI', 'Subject'])

    return c22Data, labelColumn, newSubjCount

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def addMultiIndex(c22Data, featNames, subjCount):
    ''' Adds a multi-level index to the catch 22 feature matrix data. Returns data
    as a pandas dataframe. Returns number of regions and features.'''
# NB original c22Data ndarray format: each row is a feature vector. Features from the same
# region across every subject are grouped together, i.e. First [no. of subjects]
# rows belong to one region, then the next region etc.
    [subjRows, featCols] = c22Data.shape

    roiCount = int(subjRows/subjCount)
    featCount = featCols

    # Add the multi-level index to the data
    roiIndices = list(range(0,roiCount))
    subjIndices = list(range(0,subjCount))
    multiIndex = pd.MultiIndex.from_product([roiIndices, subjIndices], names=['ROI', 'Subject'])

    c22Data = pd.DataFrame(data=c22Data, index=multiIndex, columns=featNames).rename_axis('Feature', axis=1)
    # c22Data.fillna(value=-99999,inplace=True)
    return c22Data, roiCount, featCount
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
def readLabelColumn(filePaths=[],participants=None,dataset='UCLA',indices=[]):
    ''' Returns a label vector (ndarray) with binary values where 0 = control, 1 = SCZ,
    read from either file paths or particpants.csv file, depending on dataset.'''

    if dataset == 'UCLA':
        # String '-5' is only found in all filenames that contain SCZ data
        identifier = '-5'

        labelCol = np.zeros(len(filePaths)).astype(int)
        labelCol = np.reshape(labelCol,(len(filePaths),1))

        i = 0
        for path in filePaths:
            if identifier in path:
                labelCol[i] = 1
            i += 1

    if dataset == 'COBRE':
        labelCol = participants['Diagnosis'][indices]
        labelCol.replace(to_replace=1,value=0,inplace=True) #in participants.csv, 1 indicated control
        labelCol.replace(to_replace=2,value=1,inplace=True) #in participants.csv, 2 indicated Schiz
        labelCol = np.asarray(labelCol,dtype=np.int).reshape(len(labelCol),1)

    return labelCol
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
def groupCounts(labelColumn):
    ''' Returns number of control subjects & SCZ subjects, the total number of
    subjects and the ratio of SCZ:control subjects inside a label vector. '''

    ctrlCount = (labelColumn == 0).sum()
    sczCount = (labelColumn == 1).sum()
    total = len(labelColumn)
    sczToCtrlRatio = '{0:.2f}'.format(sczCount/ctrlCount)
    return ctrlCount, sczCount, total, sczToCtrlRatio
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
def giveMeLowestThreshFD(filePath, subPath, dataset, k, threshold_fdArray):
    '''OLD FUNCTION - NOT USED'''
    ''' This function calculates the lowest threshold_fd that can be applied to still have meaningful k-folds '''

    # Initialise an index
    i = 0

    for threshold_fd in threshold_fdArray:

        # Read in and store the framewise displacement (fd) for the given dataset in a variable called fdAvgs,
        # and create the TS_path_names and indices2Keep variables

        # Store the fdAvgs
        fdAvgs = pd.read_csv(filePath,header=None);

        # Need to alphabetise and store the subject file names into a variable
        TS_path_names = sorted(glob.glob(subPath + '*.mat'))

        # Filter the subjects based on their fd, and retain the subjects that have an fd < threshold_fd
        TS_path_names, indices2Keep = removePathNames(filePath, threshold_fd, TS_path_names)
        indices2Keep = indices2Keep.tolist()

        # Adding 1 to every element in the array to convert to MATLAB indexing
        indices2KeepMat = list(np.asarray(indices2Keep) + 1)

        # Create the target column - unique for each dataset
        if dataset == 'UCLA':

            # Creating the target column
            targetCol = getTargetCol(TS_path_names)

        elif dataset == 'COBRE':

            # Creating the target column
            csvPath = '/Users/preethompal/Dropbox (Sydney Uni Student)/COBRE/participants.csv'
            COBRE = pd.read_csv(csvPath,header=None);

            targetCol = COBRE.iloc[1:,2]
            targetCol = targetCol.tolist()
            targetCol = pd.DataFrame(data=targetCol, columns=['target'])

            targetCol = targetCol.iloc[indices2Keep,:]
            targetCol = np.asarray(targetCol,dtype=np.int)

            # A '0' indicates a control subject and a '1' indicates a subject with SCZ
            targetColModified = np.where(targetCol==1, 0, targetCol) # First change the pre-existing 1s to 0s
            targetCol = np.where(targetCol==2, 1, targetColModified) # Then change the 2s to 1s

        # Store and print the subject numbers within the dataset
        Control, SCZ, Total, SCZ2Ctrl = giveMeSubjectNums(targetCol)

        if SCZ > k:
            # print(i)
            i = i+1
        elif SCZ <= k:
            lowestThreshFD = threshold_fdArray[i-1]
            print('Dataset = ' + str(dataset))
            print('Lowest Threshold FD for ' + str(k) + '-fold CV = ' + str('{0:.2f}'.format(lowestThreshFD)))
            return
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
def featureSlice(roiCount, subjCount, c22Data, featureName, subjIndicesBelowThresh):
    ''' This function takes the c22Data, feature name and the indices to be kept
    (based on the threshold fd) and returns the selected featSlice as a dataframe '''

    roiIndices = list(range(0,roiCount))
    subjIndices = list(range(0,subjCount))

    featSlice = c22Data.loc[:,featureName].values.reshape(roiCount,subjCount).transpose() # feature slice (ndarray); subjects in rows, regions in columns.

    featSlice = pd.DataFrame(data=featSlice, index=subjIndices, columns=roiIndices).rename_axis('ROI', axis=1)
    featSlice.index.name = 'Subject'
    # Return slice only with specified subjects below threshold
    return featSlice.loc[subjIndicesBelowThresh,:]
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
def tenFoldCVScore(X,y):
    ''' This function returns a 10-fold CV score after balancing the classes '''

    # Support vector machine classifier.
    svclassifier = SVC(kernel='linear', class_weight='balanced')

    # Split the data into training and test sets
    skf = StratifiedKFold(n_splits=10)
    rskf = RepeatedStratifiedKFold(n_splits=5,n_repeats=100,random_state=42)

    scores = np.zeros(500)
    i = 0
    for trainIndices, testIndices in rskf.split(X,y):
        trainIndices = trainIndices.tolist()
        testIndices = testIndices.tolist()

        X_train, X_test = X.iloc[trainIndices,:], X.iloc[testIndices,:]
        y_train, y_test = y[trainIndices], y[testIndices]
        try:
            svclassifier.fit(X_train, y_train)
        except ValueError as e:
            print("ValueError while doing tenfold CV and learn. Showing classes in y.")
            print(np.unique(y_train))
            print("Error message:", e)
            exit()
        predictions = svclassifier.predict(X_test)
        accuracy = balanced_accuracy_score(y_test, predictions)
        scores[i] = '{0:.2f}'.format(accuracy*100)
        i += 1
    return scores
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
def giveMeSVMWeights(X,y):
    ''' This function returns the feature weights when given X and y '''

    # Import the support vector classifier and balance the classes
    from sklearn.svm import SVC
    svclassifier = SVC(kernel='linear', class_weight='balanced')#, C=1e-2)

    svclassifier.fit(X, y)

    svmWeights = svclassifier.coef_
    # svmW_shape = svmWeights.shape

    return svmWeights[0]
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
def tAndPValues(labelColumn, slice):
    '''slice is either a feature or a region.'''
    '''Returns in a data frame two values (t-value from two tailed t-test, then
    p-value) in each row, one row per each ROI or feature. Also returns slice
    indices that have the most significant p-values.'''
    '''These t statistics represent (when looking at a feature for example) the
    deviation from the null hypothesis that the average of values for a feature
    in each region will be the same for schizophrenics and controls.  '''

    [subjRows, selectCols] = slice.shape
    #NB selectCols would be each roiCols if looking at a feature slice, or featCols if looking at a region slice

    tAndPVals = np.zeros([selectCols, 2])
    ctrlIndices = np.where(labelColumn == 0)[0]
    sczIndices = np.where(labelColumn == 1)[0]

    for selectCol in range(selectCols):

        ctrls = slice.iloc[ctrlIndices,selectCol]
        sczs = slice.iloc[sczIndices,selectCol]
        tAndPVal = stats.ttest_ind(ctrls, sczs, equal_var=False)
        if math.isnan(float(tAndPVal[0])):
            # print('nan made in ' + str(selectCol))
            # print(np.var(ctrls)==0)
            # print(np.var(sczs)==0)
            continue
        for j in range(2):
            tAndPVals[selectCol,j] = tAndPVal[j]

    # Since it is a two-tailed t-test, need to multiply the p-values by two (second column) [????]
    # tAndPVals[:,1] = tAndPVals[:,1] * 2

    tpValsDf = pd.DataFrame(data=tAndPVals, columns=['t-value', 'p-value'])
    tpValsDf.index.name = 'Feature / ROI'

    tpValsDf_sorted = tpValsDf.abs().sort_values(by='p-value',ascending=True)
        # features or regions with lowest p-values at the top
    indexVals = tpValsDf_sorted.index.values
    bestPValInds = indexVals[:5] #indices to the best 5 features or regions.
    return tpValsDf, bestPValInds
#-------------------------------------------------------------------------------
def tValueHistogramsOverlay(dataframesSuperList, featCount):
    fig = plt.figure()
    fig.suptitle('DiCER UCLA', y = 0.93)
    colours = sns.color_palette("cubehelix", len(dataframesSuperList))
    colours = ['r','b','g']
    for feature in range(featCount):
        ax_curr = fig.add_subplot(4,6,feature+1)
        ax_curr.set_title('Feat. ' + str(feature+1), fontsize=10)
        ax_curr.set_xlabel('t-stat')
        ax_curr.set_ylabel('No. of Regions')
        # ax_curr.legend(['procMeth1', 'procMeth2', 'procMeth3'])
        i=0
        for tpValsDfs in dataframesSuperList:
            feat_tvals = tpValsDfs[feature]['t-value']
            sns.distplot(feat_tvals, label=('procMeth'+str(i+1)), kde=False, ax=ax_curr, color=colours[i],
                        kde_kws={'color':colours[i], 'lw':2},
                        hist_kws={'histtype': 'step', 'linewidth':3})
            i+=1

    fig.subplots_adjust(left=0.07,bottom=0.07,right=0.97,top=0.88,hspace=0.55, wspace=0.33)
    handles, labels = ax_curr.get_legend_handles_labels()
    fig.legend(handles, labels)
    plt.show()


#-------------------------------------------------------------------------------
def tValueHistograms(dataset, featCount, featNames, roiCount, subjCount, c22Data, subjIndicesBelowThresh, labelColumn, returnDFs=False):
    ''' This function plots the t-distribution for each of the feature slices '''
    if not returnDFs:
        fig = plt.figure()
        fig.suptitle(str(dataset), y=0.93)
    tpValDfs = []
    for feature in range(featCount):
        featureName = featNames[feature]

        featSlice = featureSlice(roiCount,subjCount,c22Data,featureName,subjIndicesBelowThresh)
        tpValDf = tAndPValues(labelColumn,featSlice)[0]
        if returnDFs:
            tpValDfs.append(tpValDf)
            continue

        ax_curr= fig.add_subplot(4,6,feature+1)

        feat_tvals = tpValDf['t-value']

        # Plot the histograms
        sns.distplot(feat_tvals, kde=False, fit=stats.norm, ax=ax_curr)
        ax_curr.set_title('Feat. ' + str(feature+1), fontsize=10)
        ax_curr.set_xlabel('t-stat')
        ax_curr.set_ylabel('Frequency')

    if returnDFs:
        return tpValDfs

    fig.subplots_adjust(left=0.07,bottom=0.07,right=0.97,top=0.88,hspace=0.55, wspace=0.33)
    plt.show()
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
def plotPCA(slice, labelColumn):
    ''' This function displays a figure derived from PCA '''
    # Standardizing the features
    X = StandardScaler().fit_transform(slice)

    pca = PCA(n_components=2)
    principalComponents = pca.fit_transform(X)
    principalDf = pd.DataFrame(data=principalComponents, columns=['PC1', 'PC2'])

    targetCol_df = pd.DataFrame(data=labelColumn, columns=['Diagnosis'])

    finalDf = pd.concat([principalDf, targetCol_df], axis = 1)

    # Plotting the PCA scatterplot using Seaborn
    sns.set()
    ax = sns.relplot(x='PC1',y='PC2',data=finalDf,hue='Diagnosis',palette='Set1')
    plt.show()
    return
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
def displayViolinPlots(labelColumn, bestPValInds, slice, selection, selectionNum):
    '''Displays violin plots of five feature columns specified by bestPValInds - the
    learning features with the lowest p-values from two-tailed t test that compared
    feature columns across schiz and controls.'''

    subplotIndex = 1;
    fig = plt.figure()

    controlRows = np.where(labelColumn == 0)[0]
    schizRows = np.where(labelColumn == 1)[0]

    # Z-score the data
    sliceDF = pd.DataFrame(data=slice)
    sliceDF = sliceDF.apply(zscore)

    for i in bestPValInds:
        goodFeatColumnCtrls = sliceDF.iloc[controlRows,i]
        goodFeatColumnSczs = sliceDF.iloc[schizRows,i]

        goodFeatDF = pd.DataFrame({'SCZ':goodFeatColumnSczs,'Control':goodFeatColumnCtrls})

        # Violin plots
        ax = fig.add_subplot(2,3,subplotIndex)
        ax = sns.violinplot(data=goodFeatDF, order=['SCZ','Control'])
        plt.xlabel('Diagnosis')

        if selection == 0:
            ylabel = 'ROI ' + str(i)
        elif selection == 1:
            ylabel = 'Feature ' + str(i)
        plt.ylabel(ylabel)

        subplotIndex += 1;

    if selection == 0:
         plt.suptitle('Top 5 ROI in Feature ' + str(selectionNum), fontsize=12.5, y=1.05)
    elif selection == 1:
        plt.suptitle('Top 5 Features in ROI ' + str(selectionNum), fontsize=12.5, y=1.05)
    plt.tight_layout()
    plt.show()
    return
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def accuracyOfRegionsAndFeatures(c22Data, labelColumn, roiCount, subjCount, featCount):
    regionsAndFeaturesMatrix = np.zeros((subjCount, featCount*roiCount))
    for subj in range(subjCount):
        for roi in range(roiCount):
            regionsAndFeaturesMatrix[subj,roi*featCount:(roi+1)*featCount] = c22Data.loc[roi, subj]

    regionsAndFeaturesMatrix = pd.DataFrame(data=regionsAndFeaturesMatrix, index=range(subjCount))
    regionsAndFeaturesMatrix.rename_axis('Subjects',inplace=True)
    regionsAndFeaturesMatrix.rename_axis('RegionsFeatures', axis='columns',inplace=True)
    print(regionsAndFeaturesMatrix)
    # not applying zcore.
    # features x regions. subjects stacked on top of eachother.

    regionsAndFeaturesMatrixZScored = regionsAndFeaturesMatrix.apply(zscore)

    for featureRegion, featCol in regionsAndFeaturesMatrixZScored.iteritems():
        if featCol.isnull().values.any():
            print('Removing featureRegion '+str(featureRegion)+" for cross validation due to presence of NaNs in this feature's z-scores.")
            print("Presence of NaNs: ", str(np.mean(featCol.isnull().values)*100)+'%')
            regionsAndFeaturesMatrixZScored.drop(featureRegion, axis='columns', inplace=True)

    X = regionsAndFeaturesMatrixZScored
    y = np.ravel(labelColumn)

    tenFoldScore = tenFoldCVScore(X, y)
    accuracy = tenFoldScore.mean()
    error = tenFoldScore.std()

    print('Mean Accuracy (across all regions and features at once) = ' +
    "{0:.2f}".format(accuracy) + '% +/- ' + "{0:.2f}".format(error))

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def accuracyOfRegions(roiCount, c22Data, subjIndicesBelowThresh, labelColumn, dispFigs, returnDF=False):
    ''' By default returns the mean accuracy and error of all regions. Turning
    on dispFigs will show the distribution of accuracies for each region, and
    which five regions were the most accurate. Can optionally return the data
    frame of all region accuracies/errors.'''

    roiNums = np.zeros([roiCount])
    roiAccuracies = np.zeros([roiCount])
    roiErrors = np.zeros([roiCount])

    for n in range(0, roiCount):

        c22DataRegionSlice = c22Data.loc[n,subjIndicesBelowThresh,:]
        c22DataRegionSliceZScored = c22DataRegionSlice.apply(zscore)

        for featName, featCol in c22DataRegionSliceZScored.iteritems():
            if featCol.isnull().values.any():
                print('Removing feature '+featName+' for cross validation during region '+str(n) +" due to presence of NaNs in this feature's z-scores.")
                print("Presence of NaNs: ", str(np.mean(featCol.isnull().values)*100)+'%')
                c22DataRegionSliceZScored.drop(featName, axis='columns', inplace=True)
        X = c22DataRegionSliceZScored
        y = np.ravel(labelColumn)

        # Denotes how good all the features in a given region are at predicting if a subject has SCZ
        sliceTenFoldScore = tenFoldCVScore(X,y)
        roiNums[n] = n
        roiAccuracies[n] = sliceTenFoldScore.mean()
        roiErrors[n] = sliceTenFoldScore.std()

    df = pd.DataFrame({'Region':roiNums,'% Accuracy':roiAccuracies,'ROI Error':roiErrors})
    df['Region'] = df.Region.astype(int)
    df.set_index('Region', inplace=True)
    # df.drop(problemRegions, inplace=True)

    df_sorted = df.sort_values(by='% Accuracy',ascending=False)

    if dispFigs:
        print('\n'+df_sorted[:5].to_string()+'\n')

        print('Mean Accuracy (across all regions) = ' +
        "{0:.2f}".format(df['% Accuracy'].mean()) + '%\n')

        print('Mean Error (across all regions) = ' +
        "{0:.2f}".format(df['ROI Error'].mean()) + '%\n')
        # roiAccuracies = np.delete(np.asarray(roiAccuracies), problemRegions)
        # sns.distplot(roiAccuracies, bins='auto', kde=False)
        # plt.xlabel('Classification Accuracy (%)')
        # plt.ylabel('Number of Regions')
        # plt.show()

    if returnDF:
        return df
    meanROIAcc = '{0:.2f}'.format(df['% Accuracy'].mean())
    meanROIError = '{0:.2f}'.format(df['ROI Error'].mean())
    return meanROIAcc, meanROIError
#-------------------------------------------------------------------------------
def distributionsOverlayPlot(dataframes, columnName, plotName):
    colours = sns.color_palette("cubehelix", len(dataframes))
    fig, ax = plt.subplots()
    fig.suptitle(plotName)
    i = 0
    for df in dataframes:
        plotValues = np.asarray(df[columnName])
        print('\nMean (across all '+plotName+') = ' +
        "{0:.2f}".format(df[columnName].mean()) + '%\n')
        sns.distplot(plotValues, bins='auto', kde=False, color=colours[i], label='procMeth'+str(i+1), hist_kws={'histtype': 'step', 'linewidth':3})
        i+=1
    plt.xlabel(plotName + ' Classification Accuracy (%)')
    plt.ylabel('Number of '+plotName)
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels)
    plt.show()

#-------------------------------------------------------------------------------
def accuracyOfFeatures(c22Data, roiCount, subjCount, featNames, subjIndicesBelowThresh, labelColumn, dispFigs, featCount=22, returnDF=False):
    ''' By default returns the mean accuracy and error of all features. Turning
    on dispFigs will show the distribution of accuracies for each feature, and
    which five features were the most accurate. Can optionally return the sorted
    (by accuracy) data frame of all feature accuracies/errors.'''

    featNums = np.zeros([featCount])
    featAccuracies = np.zeros([featCount])
    featErrors = np.zeros([featCount])

    for n in range(0, featCount):
        # For each of the 1 to 22 features, show me how good all of the ROIs are at predicting
        # whether the subject has SCZ or not
        c22DataFeatureSlice = featureSlice(roiCount, subjCount, c22Data, featNames[n], subjIndicesBelowThresh)
        c22DataFeatureSliceZScored = c22DataFeatureSlice.apply(zscore)

        for regNum, regCol in c22DataFeatureSliceZScored.iteritems():
            if regCol.isnull().values.any():
                print('Removing region '+str(regNum)+' for cross validation during feature '+str(n) +" due to presence of NaNs in this region's z-scores.")
                print("Presence of NaNs: ", str(np.mean(regCol.isnull().values)*100)+'%')
                c22DataFeatureSliceZScored.drop(regNum, inplace=True, axis='columns')

        # Classification: features and labels.
        X = c22DataFeatureSliceZScored
        y = np.ravel(labelColumn)

        sliceTenFoldScore = tenFoldCVScore(X,y)
        featNums[n] = n
        featAccuracies[n] = sliceTenFoldScore.mean()
        featErrors[n] = sliceTenFoldScore.std()

    df = pd.DataFrame({'Feature':featNums,'% Accuracy':featAccuracies,'Feat Error':featErrors})
    df['Feature'] = df.Feature.astype(int)
    df.set_index('Feature',inplace=True)
    df_sorted = df.sort_values(by='% Accuracy',ascending=False)


    if dispFigs:
        print(df_sorted[:5].to_string())

        print('\nMean Accuracy (across all features) = ' +
        "{0:.2f}".format(df['% Accuracy'].mean()) + '%\n')

        print('Mean Error (across all features) = ' +
        "{0:.2f}".format(df['Feat Error'].mean()) + '%\n')
        sns.distplot(featAccuracies, bins='auto', kde=False)
        plt.xlabel('Classification Accuracy (%)')
        plt.ylabel('Number of Features')
        plt.show()

    meanFeatAcc = '{0:.2f}'.format(df['% Accuracy'].mean())
    meanFeatError = '{0:.2f}'.format(df['Feat Error'].mean())
    if returnDF:
        return df
    return meanFeatAcc, meanFeatError
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
def regByRegAnalysis(roi, c22Data, labelColumn, roiCount, subjIndicesBelowThresh):
    ''' Shows accuracy of classification when learning from a single specified region. '''

    roiSlice = c22Data.loc[roi,subjIndicesBelowThresh,:]
    roiSliceZScored = roiSlice.apply(zscore)

    X = roiSliceZScored
    y = np.ravel(labelColumn)

    sliceTenFoldScore = tenFoldCVScore(X,y)
    avgScore = '{0:.2f}'.format(sliceTenFoldScore.mean())
    avgSTD = '{0:.2f}'.format(sliceTenFoldScore.std())

    print('Analysis of Region ' + str(roi) + ':\n')
    print("Accuracy as a percentage: " + avgScore + " +/- " + avgSTD)

    #Use tAndPValues to get indices to 5 features with the most significant p values.
    tpValsDf, bestPValInds = tAndPValues(labelColumn, roiSlice)

    # Show me the PCA figure
    plotPCA(roiSlice, labelColumn)

    # Show me the top five features as violin plots in the ROI being analysed
    displayViolinPlots(labelColumn, bestPValInds, roiSlice, 1, roi)

    return avgScore, avgSTD
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
def featByFeatAnalysis(featureName, c22Data, labelColumn, roiCount, subjCount, subjIndicesBelowThresh):
    ''' Shows accuracy of classification when learning from a single specified feature. '''

    featSlice = featureSlice(roiCount, subjCount, c22Data, featureName, subjIndicesBelowThresh)
    featSliceZScored = featSlice.apply(zscore)

    for regNum, regCol in featSliceZScored.iteritems():
        if regCol.isnull().values.any():
            print('Removing region '+str(regNum)+' for cross validation during feature '+featureName+" due to presence of NaNs in this region's z-scores.")
            print("Presence of NaNs: ", str(np.mean(regCol.isnull().values)*100)+'%')
            featSliceZScored.drop(regNum, inplace=True, axis='columns')

    X = featSliceZScored
    y = np.ravel(labelColumn)

    sliceTenFoldScore = tenFoldCVScore(X,y)
    avgScore = '{0:.2f}'.format(sliceTenFoldScore.mean())
    avgSTD = '{0:.2f}'.format(sliceTenFoldScore.std())

    print('Analysis of Feature ' + str(featureName) + ':\n')
    # print('10-fold CV scores as a percentage: ' + str(scores)+'\n')

    # Mean 10-fold CV score with an error of 1 std dev
    print("Accuracy as a percentage: " + avgScore + " +/- " + avgSTD)

    #Use tAndPValues to get indices to 5 regions with the most significant p values.
    tpValsDf, bestPValInds = tAndPValues(labelColumn, featSlice)

    # Show me the PCA figure
    plotPCA(featSlice, labelColumn)

    # Show me the top five ROIs as violin plots in the feature being analysed
    displayViolinPlots(labelColumn, bestPValInds, featSlice, 0, featureName)

    return avgScore, avgSTD
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
def featurePerformanceJointPlot(sortedFeats,sortedFeatsAlt, altName):
    ''' This function plots the jointplot of the two datasets being compared,
    namely UCLA and COBRE '''
    if altName == 'UCLA':
        df1Name = 'COBRE'
    if altName == 'COBRE':
        df1Name = 'UCLA'

    df1 = sortedFeats;
    df1 = df1.sort_values(by='Feature')
    df1.columns = [df1Name+' Feat Acc (%)', df1Name+' Feat Error (%)']
    df1FeatAccs = df1.iloc[:,0]
    df1SDs = df1.iloc[:,1]
    df1 = df1.rename_axis(df1Name+' Feature')

    print(df1)
    print('')

    df2 = sortedFeatsAlt;
    df2 = df2.sort_values(by='Feature')
    df2.columns = [altName+' Feat Acc (%)', altName+' Feat Error (%)']
    df2FeatAccs = df2.iloc[:,0]
    df2SDs = df2.iloc[:,1]
    df2 = df2.rename_axis(altName+' Feature')

    print(df2)

    sns.jointplot(df1FeatAccs, df2FeatAccs, kind='reg')
    plt.show()
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
def jointAccNullDistributionPlot(accuracies, dataframe):

    pass
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------


def kiloLabelShufflesAndLearnsRegionsJoint(labelColumn, c22Data, roiCount):
    df = pd.DataFrame({'Iteration':[], 'Average Accuracy': []})
    for i in range(1000):
        individualAccuracies = np.zeros(roiCount)
        for roi in range(roiCount):
            roiSlice = c22Data.loc[roi,:,:]
            roiSliceZScored = roiSlice.apply(zscore)
            for featName, featCol in roiSliceZScored.iteritems():
                if featCol.isnull().values.any():
                    print('Removing feature '+featName+' for cross validation during region '+str(roi) +" due to presence of NaNs in this feature's z-scores.")
                    print("Presence of NaNs: ", str(np.mean(featCol.isnull().values)*100)+'%')
                    roiSliceZScored.drop(featName, axis='columns', inplace=True)

            X = roiSliceZScored
            np.random.shuffle(labelColumn)
            y = np.ravel(labelColumn)
            sliceTenFoldScore = tenFoldCVScore(X,y)
            meanScore = sliceTenFoldScore.mean()
            individualAccuracies[roi] = meanScore
        df = df.append({'Iteration':i,'Average Accuracy': individualAccuracies.mean()},ignore_index=True)

    return df
#-------------------------------------------------------------------------------
def kiloLabelShufflesAndLearnsRegions(labelColumn, c22Data, roiCount):
    df = pd.DataFrame({'Region':[], 'Iteration': [], 'Average Accuracy': []})

    for roi in range(roiCount):
        if roi == 0:
            continue
        roiSlice = c22Data.loc[roi,:,:]
        roiSliceZScored = roiSlice.apply(zscore)
        print("kiloLabelShufflesAndLearnsRegions is doing roi", roi)

        for featName, featCol in roiSliceZScored.iteritems():
            if featCol.isnull().values.any():
                print('Removing feature '+featName+' for cross validation during region '+str(roi) +" due to presence of NaNs in this feature's z-scores.")
                print("Presence of NaNs: ", str(np.mean(featCol.isnull().values)*100)+'%')
                roiSliceZScored.drop(featName, axis='columns', inplace=True)
        X = roiSliceZScored
        for i in range(275):
            np.random.shuffle(labelColumn)
            y = np.ravel(labelColumn)
            sliceTenFoldScore = tenFoldCVScore(X,y)
            meanScore = sliceTenFoldScore.mean()
            df = df.append({'Region':roi, 'Iteration':i,'Average Accuracy':meanScore},ignore_index=True)

    return df
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
def kiloLabelShufflesAndLearnsFeatures(labelColumn, c22Data, subjIndicesBelowThresh, roiCount, subjCount, featNames, featCount=22):

    df = pd.DataFrame({'Feature': [], 'Iteration': [], 'Average Accuracy': [], 'stdDev': []})

    for feature in range(featCount):
        featName = featNames[feature]
        featSlice = featureSlice(roiCount, subjCount, c22Data, featName, subjIndicesBelowThresh)
        featSliceZScored = featSlice.apply(zscore)
        for regNum, regCol in featSliceZScored.iteritems():
            if regCol.isnull().values.any():
                print('Removing region '+str(regNum)+' for cross validation during feature '+str(feature) +" due to presence of NaNs in this region's z-scores.")
                print("Presence of NaNs: ", str(np.mean(regCol.isnull().values)*100)+'%')
                featSliceZScored.drop(regNum, inplace=True, axis='columns')

        X = featSliceZScored
        for i in range(1000):
            np.random.shuffle(labelColumn)
            y = np.ravel(labelColumn)
            sliceTenFoldScore = tenFoldCVScore(X,y)
            meanScore = sliceTenFoldScore.mean()
            scoreStdDev = sliceTenFoldScore.std()
            df = df.append({'Feature':feature, 'Iteration':i,'Average Accuracy':meanScore,'stdDev':scoreStdDev},ignore_index=True)

    return df
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def jointROIAccuracyPValTriple(accuracies, randomLearnData, roiCount):
    pass
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def roiAccuracyPValsTriple(randomLearnData, roiAccsDataframes, roiCount):
    pVals = [np.zeros(roiCount) for i in range(3)]

    roiAccs = []
    for df in roiAccsDataframes:
        roiAccs.append(np.asarray(df['% Accuracy']))

    if randomLearnData.index.name != 'Region':
        randomLearnData.set_index('Region', inplace=True)

    randomAccsAll = np.asarray(randomLearnData['Average Accuracy'])

    for roi in range(roiCount):
        for i in range(3):
            pVals[i][roi] = np.mean(randomAccsAll>=roiAccs[i][roi])

    pValsCorrected = []
    for i in range(3):
        pValsCorrected.append(sm.stats.multitest.multipletests(pVals[i], method='fdr_bh')[1])

    df = pd.DataFrame({'procMeth1':pVals[0], 'procMeth1 Corrected':pValsCorrected[0],
                        'procMeth2':pVals[1], 'procMeth2 Corrected':pValsCorrected[1],
                        'procMeth3':pVals[2], 'procMeth3 Corrected':pValsCorrected[2]})
    df = df.round(3)
    df.index.name = 'Region'
    df.index += 1
    pd.options.display.width = 0
    pd.set_option('display.max_rows', None)
    print(df)

    for index, row in df.iterrows():
        significant = False
        for pVal in list([row[1],row[3],row[5]]):
            if pVal <= 0.05:
                significant = True
        if not significant:
            df.drop(index, inplace=True)

    print(df)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def jointFeatureAccuracyPValTriple(accuracies, randomLearnData, featCount):
    pass
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def featureAccuracyPValsTriple(featAccDataframes, randomLearnData, featCount):
    ''' Prints the combined p-values (probability of randomly calculated accuracy
     being equal to or greater than correctly calculated accuracy) for both datasets '''
    pVals = [np.zeros(featCount) for i in range(3)]
    featAccs = []
    for df in featAccDataframes:
        featAccs.append(df['% Accuracy'])

    if randomLearnData.index.name != 'Feature':
        randomLearnData.set_index('Feature',inplace=True)

    randomAccsAll = np.asarray(randomLearnData['Average Accuracy'])

    for feature in range(featCount):

        # What is the likelihood that a randomly calculated classification
        # accuracy will beat the mean feature accuracy that has been calculated?
        for i in range(3):
            pVals[i][feature] = np.mean(randomAccsAll>=featAccs[i][feature])

    # Multiple hypothesis test correction. Benjamini/Hochberg (non-negative)
    pValsCorrected = []
    for i in range(3):
        pValsCorrected.append(sm.stats.multitest.multipletests(pVals[i], method='fdr_bh')[1])

    df = pd.DataFrame({'procMeth1':pVals[0], 'procMeth1 Corrected':pValsCorrected[0],
                        'procMeth2':pVals[1], 'procMeth2 Corrected':pValsCorrected[1],
                        'procMeth3':pVals[2], 'procMeth3 Corrected':pValsCorrected[2]})
    df = df.round(3)
    df.index.name = 'Feature'
    pd.options.display.width = 0
    print(df)
    # Combining...
    # Note: not using corrected raw p values for each dataset before combining. Correcting again after combining.
    pValueTrips = np.column_stack((pVals[0], pVals[1], pVals[2]))
    pValsCombined = []
    for pVals in pValueTrips:
        pValsCombined.append(stats.combine_pvalues(pVals, method='fisher')[1])

    # Multiple hypothesis test correction. Benjamini/Hochberg (non-negative)
    pValsCombinedCorrected = sm.stats.multitest.multipletests(pValsCombined, method='fdr_bh')[1]
    df = pd.DataFrame({'Combined': pValsCombined, 'Corrected': pValsCombinedCorrected})
    df = df.round(3)
    df.index.name = 'Feature'
    print(df)
#-------------------------------------------------------------------------------
def featureAccuracyPVals(dfFeatAccsUCLA,dfFeatAccsCOBRE,randomLearnDataUCLA,randomLearnDataCOBRE,featCount):
    ''' Prints the combined p-values (probability of randomly calculated accuracy
     being equal to or greater than correctly calculated accuracy) for both datasets '''
    print('UCLA',dfFeatAccsUCLA)
    # print('COBRE',dfFeatAccsCOBRE)

    pValsUCLA = np.zeros(featCount)
    # pValsCOBRE = np.zeros(featCount)

    featAccsUCLA = np.asarray(dfFeatAccsUCLA.iloc[:,0])
    uclaSDs = np.asarray(dfFeatAccsUCLA.iloc[:,1])

    # featAccsCOBRE = np.asarray(dfFeatAccsCOBRE.iloc[:,0])
    # cobreSDs = np.asarray(dfFeatAccsCOBRE.iloc[:,1])

    randomLearnDataUCLA.set_index('Feature',inplace=True)
    # randomLearnDataCOBRE.set_index('Feature',inplace=True)

    for feature in range(featCount):
        if feature==3:
            continue

        randomLearnSliceUCLA = randomLearnDataUCLA.loc[feature]
        # randomLearnSliceCOBRE = randomLearnDataCOBRE.loc[feature]

        randomAccsUCLA = np.asarray(randomLearnSliceUCLA['Average Accuracy'])
        # randomAccsCOBRE = np.asarray(randomLearnSliceCOBRE['Average Accuracy'])

        # What is the likelihood that a randomly calculated classification
        # accuracy will beat the mean feature accuracy that has been calculated?

        pValsUCLA[feature] = np.mean(randomAccsUCLA>=featAccsUCLA[feature])
        # pValsCOBRE[feature] = np.mean(randomAccsCOBRE>=featAccsCOBRE[feature])

    # Multiple hypothesis test correction. Benjamini/Hochberg (non-negative)
    pValsUCLA = np.delete(pValsUCLA, 3)
    pValsCorrectedUCLA = sm.stats.multitest.multipletests(pValsUCLA, method='fdr_bh')[1]
    # pValsCorrectedCOBRE = sm.stats.multitest.multipletests(pValsCOBRE, method='fdr_bh')[1]
    df = pd.DataFrame({'UCLA P Values':pValsCorrectedUCLA})
    df = df.round(3)
    idx = np.asarray([i for i in range(22)])
    idx = np.delete(idx,3)
    df.index = idx
    df.index.name = 'Feature'
    print(df)
    # Note: not using corrected raw p values for each dataset before combining. Correcting again after combining.
    # pValuePairs = np.column_stack((pValsUCLA, pValsCOBRE))
    # pValsCombined = []
    # for pVals in pValuePairs:
    #     pValsCombined.append(stats.combine_pvalues(pVals, method='fisher')[1])
    # # Multiple hypothesis test correction. Benjamini/Hochberg (non-negative)
    # pValsCombined = sm.stats.multitest.multipletests(pValsCombined, method='fdr_bh')[1]
    # df = pd.DataFrame(data=pValsCombined, columns=['Combined P-Values'])
    # df = df.round(3)
    # df.index.name = 'Feature'
    # print(df)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def roiAccNullDistributionsPlotTriple(randomLearnData, roiAccsDataframes, roiCount):
    roiAccs = []
    for df in roiAccsDataframes:
        roiAccs.append(np.asarray(df['% Accuracy']))

    randomLearnData.set_index('Region',inplace=True)

    randomAccsAll = randomLearnData['Average Accuracy']
    colours = ['r','b','g']

    i = 1
    figure = plt.figure(i)
    figure.suptitle('Individual Region Accuracy vs Null Distribution')
    adjusted = False
    pos = 1
    for roi in range(roiCount):
        if (roi%24==0) and roi != 0:
            pos = 1
            i+=1
            figure = plt.figure(i)
            figure.suptitle('Individual Region Accuracy vs Null Distribution')
            adjusted = False
        ax = figure.add_subplot(4,6, pos)
        pos+=1
        ax.set_xlabel('Accuracy')
        ax.set_ylabel('Frequency')
        ax.set_title('Region '+str(roi+1))
        ax.hist(randomAccsAll,bins=50)
        for j in range(3):
            ax.axvline(x=roiAccs[j][roi], color=colours[j], linestyle='dashed', label=('procMeth'+str(j+1)))

        if not adjusted:
            figure.subplots_adjust(left=0.07,bottom=0.07,right=0.97,top=0.88,hspace=0.55, wspace=0.33)
            handles, labels = ax.get_legend_handles_labels()
            figure.legend(handles, labels)
            adjusted = True

    plt.show()

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def featAccNullDistributionsPlotTriple(randomLearnData, featAccsDataframes, featCount):

    featAccs = []
    for df in featAccsDataframes:
        featAccs.append(np.asarray(df['% Accuracy']))

    randomLearnData.set_index('Feature',inplace=True)

    figure = plt.figure(1)
    figure.suptitle('UCLA DiCER')

    randomAccsAll = randomLearnData['Average Accuracy']
    colours = ['r','b','g']

    for feature in range(featCount):
        ax = figure.add_subplot(4,6,feature+1)
        ax.set_xlabel('Accuracy')
        ax.set_ylabel('Frequency')
        ax.set_title('Feature '+str(feature+1))
        ax.hist(randomAccsAll,bins=50)
        for i in range(3):
            ax.axvline(x=featAccs[i][feature], color=colours[i],linestyle='dashed',label=('procMeth'+str(i+1)))

    figure.subplots_adjust(left=0.07,bottom=0.07,right=0.97,top=0.88,hspace=0.55, wspace=0.33)
    handles, labels = ax.get_legend_handles_labels()
    figure.legend(handles, labels)
    plt.show()
    return
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
def featAccNullDistributionsPlot(randomLearnDataUCLA, randomLearnDataCOBRE, dfFeatAccsUCLA, dfFeatAccsCOBRE, featCount):
    featAccsUCLA = np.asarray(dfFeatAccsUCLA.iloc[:,0])
    featAccsCOBRE = np.asarray(dfFeatAccsCOBRE.iloc[:,0])

    randomLearnDataUCLA.set_index('Feature',inplace=True)
    randomLearnDataCOBRE.set_index('Feature',inplace=True)
    for i in range(2):
        figure = plt.figure(i+1)
        if i:
            figure.suptitle('COBRE')
        else:
            figure.suptitle('UCLA')

        for feature in range(featCount):
            randomAccsOneFeatureBothDatasets = [randomLearnDataUCLA.loc[feature]['Average Accuracy'],
                                                randomLearnDataCOBRE.loc[feature]['Average Accuracy']]
            ax = figure.add_subplot(4,6,feature+1)
            ax.set_xlabel('Accuracy')
            ax.set_ylabel('Frequency')
            ax.set_title('Feature '+str(feature+1))
            ax.hist(randomAccsOneFeatureBothDatasets[i],bins=50)
            if i:
                ax.axvline(x=featAccsCOBRE[feature], color='g',linestyle='dashed')
            else:
                ax.axvline(x=featAccsUCLA[feature], color='g',linestyle='dashed')
    plt.show()
    return
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
def accuraciesVaryingWithFDPlot(learnData):
    ''' This function plots the balanced classification accuracies of features
    and regions as they vary with fd'''
    avgROIError = learnData['AverageErrorUsingRegions']
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
    return
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
def showMeFDAcrossSubs(learnData,fdAvgs):
    ''' This function does two things when looking at a given dataset:
    1. It displays the fd variation across subjects
    2. It displays the number of subjects left as they are filtered by fd
    filePathB = 'fdAvgs_DATASET.txt' '''

    # fd variation across subjects
    fdArray = learnData;

    fd = fdArray.iloc[:,0]
    fdAvgs.hist(column='Avgs', bins='auto')

    plt.xlim(max(fd), min(fd))
    plt.title('FD Distribution')
    plt.xlabel('FD Averages (cm)')
    plt.ylabel('No. of Subjects')
    plt.show()

    # No. of subjects left as fd decreases

    SCZ = fdArray.iloc[:,1]
    Control = fdArray.iloc[:,2]

    fdArray.plot.line(x='fd', y=['SCZ', 'Control'])

    plt.xlim(max(fd), min(fd))
    plt.title('FD vs Subjects Remaining')
    plt.xlabel('FD Averages (cm)')
    plt.ylabel('No. of Subjects')
    plt.show()
#-------------------------------------------------------------------------------
