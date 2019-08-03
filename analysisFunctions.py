#!/usr/local/bin/python3

# Import required modules
import glob
import scipy.io as sio
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('ggplot')
import seaborn as sns
import sklearn

from scipy.stats import zscore

#-------------------------------------------------------------------------------
def removePathNames(filePath, threshold_fd, TS_path_names):
    ''' This function modifies the entries within TS_path_names based on the
    threshold fd value '''

    fdAvgs = pd.read_csv(filePath,header=None);
    indices2Del = np.where(fdAvgs > threshold_fd)[0]
    indices2Keep = np.where(fdAvgs <= threshold_fd)[0]

    for i in sorted(indices2Del, reverse=True):
        del TS_path_names[i]
    return TS_path_names, indices2Keep
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
def addIndices(element, subPath, PyFeatList):
    ''' This function adds a multi-level index to the feature matrix data that
        is read in as a text file, and returns the matrix
        The function takes two other inputs: a .txt file of the feature names
        as well as the subject file path for the given dataset '''

    from numpy import genfromtxt
    tsData = genfromtxt(element, delimiter=',')
    [rows, cols] = tsData.shape

    # Store the number of subjects and ROIs
    TS_path_names = sorted(glob.glob(subPath + '*.mat'))

    subjects = len(TS_path_names)
    ROIs = int(rows/subjects)
    feats = cols

    # Add the multi-level index to the data
    ROI_index = list(range(1,ROIs+1))
    sub_index = list(range(1,subjects+1))

    iterables = [ROI_index, sub_index]
    index = pd.MultiIndex.from_product(iterables, names=['ROI', 'Subject'])

    featList = [lines.rstrip('\n') for lines in open(PyFeatList)]
    tsData = pd.DataFrame(data=tsData, index=index, columns=featList).rename_axis('Feature', axis=1)
    return tsData, ROIs, subjects, feats, featList
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
def getTargetCol(TS_path_names):
    ''' This function, when called, prompts the user to input a string
    which can be used to uniquely identify any SCZ data file
    It returns a binary target column (0 = control, 1 = SCZ) '''

    # Enter a string found in all filenames that only contain SCZ data
    string = '-5'

    # Initialise and format a target column of length TS_path_names
    x = np.zeros(len(TS_path_names))
    x = x.astype(int)
    Diagnosis = np.reshape(x,(len(TS_path_names),1))

    # Initialise an index
    i = 0

    for path in TS_path_names:
        if string in path:
            # If the string is found, the data must be from a SCZ patient
            Diagnosis[i] = 1
            # Increment index
            i += 1
        else:
            # Increment index
            i += 1
    return Diagnosis
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
def giveMeSubjectNums(targetCol):
    ''' This function takes the target column as an input and returns the number of control subjects &
    SCZ subjects, the total number of subjects and the ratio of SCZ : control subjects (SCZ2Ctrl) '''

    Control = (targetCol == 0).sum()
    SCZ = (targetCol == 1).sum()
    Total = int(SCZ + Control)
    SCZ2Ctrl = '{0:.2f}'.format(SCZ/Control)
    return Control, SCZ, Total, SCZ2Ctrl
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
def giveMeLowestThreshFD(filePath, subPath, dataset, k, threshold_fdArray):
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
            csvPath = '/Users/AV/Dropbox/COBRE/participants.csv'
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

''' Selecting a ROI slice is achieved by the one-liner below: '''

# ROISlice = tsData.loc[ROI,indices2KeepMat,:] # The ROI needs to be chosen

#-------------------------------------------------------------------------------
def getFeatSlice(ROIs, subjects, tsData, featureName, indices2KeepMat):
    ''' This function takes the tsData, feature name and the indices to be kept
    (based on the threshold fd) and returns the selected featSlice as a dataframe '''

    ROI_index = list(range(1,ROIs+1))
    sub_index = list(range(1,subjects+1))

    featSlice = tsData.loc[:,featureName].values.reshape(ROIs,subjects).transpose() # Make the featSlice (no index)

    featSlice = pd.DataFrame(data=featSlice, index=sub_index, columns=ROI_index).rename_axis('ROI', axis=1) # Add the index
    featSlice.index.name = 'Subject' # Name the index
    featSlice = featSlice.loc[indices2KeepMat,:] # Take a subsection of the featSlice
    return featSlice
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
def get10FoldCVScore(X,y):
    ''' This function returns a 10-fold CV score after balancing the classes '''

    # Import the support vector classifier and balance the classes
    from sklearn.svm import SVC
    svclassifier = SVC(kernel='linear', class_weight='balanced')

    # Split the data into training and test sets
    from sklearn.model_selection import StratifiedKFold
    skf = StratifiedKFold(n_splits=10)

    # Import accuracy score
    from sklearn.metrics import balanced_accuracy_score

    # Initialise a few variables
    scores = np.zeros(10)
    i = 0

    for train_index, test_index in skf.split(X,y):

        train_index = train_index.tolist()
        test_index = test_index.tolist()

#         print("Train:", train_index)
#         print('')
#         print("Validation:",test_index)

        X_train, X_test = X.iloc[train_index,:], X.iloc[test_index,:]
        y_train, y_test = y[train_index], y[test_index]

        svclassifier.fit(X_train, y_train)
        y_pred = svclassifier.predict(X_test)

#         print('')
#         print('y_test = ', y_test)
#         print('')
#         print('y_pred = ', y_pred)
#         print('')

        scores[i] = '{0:.2f}'.format(balanced_accuracy_score(y_test, y_pred)*100)
#         print('Acc % = ', scores[i])
#         print('')

        # Increment index
        i += 1
    return scores
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
def getTPVals(targetCol, DataSlice):
    ''' This function computes the t-values (from a two-tail t-test) and
    the p-values
    It stores these two values (t-value, then p-value) in each row, one row per
    each ROI or feature
    It returns a few outputs, sigPValInds being the main one '''

    from scipy.stats import stats,ttest_ind

    # Initialise the array and assign its shape
    DataSlice = pd.DataFrame(data=DataSlice)
    [rows, cols] = DataSlice.shape
    tpValArray = np.zeros([cols, 2])

    # Find all instances of 0 in the target column (which would indicate control data)
    # and store the indices in an array
    a = np.where(targetCol == 0)[0]

    # Find all instances of 1 in the target column (which would indicate SCZ data)
    # and store the indices in an array
    b = np.where(targetCol == 1)[0]

    # Loop through the array, calculate and store the t and p values
    for i in range(cols):

        # Calculate the t and p values by inputting the two 'halves' of every column in the DataSlice,
        # into the ttest_ind function
        controlFeatCol = DataSlice.iloc[a,i]
        SCZFeatCol = DataSlice.iloc[b,i]

        # Store the statistics in the variable, tpVal (which changes on each iteration of the outer loop)
        tpVal = stats.ttest_ind(controlFeatCol, SCZFeatCol, equal_var=False)

        for j in range(2):

            # Store the values into each column
            tpValArray[i,j] = tpVal[j]

    # Since it is a two-tailed t-test, need to multiply the p-values by two (second column)
    tpValArray[:,1] = tpValArray[:,1] * 2

    # Formatting the tpValArray
    tpValDf = pd.DataFrame(data=tpValArray, columns=['t-value', 'p-value'])
    tpValDf.index.name = 'Feature / ROI'
    tpValDf.index += 1 # Make the starting index 1

    # Sort the data (including the indices) by p-value significance
    tpValDf_sorted = tpValDf.abs().sort_values(by='p-value',ascending=True)

    # Store the first five indices of the sorted dataframe - will need to use these
    # indices to access the relevant feature columns in X, the feature matrix
    indexVals = tpValDf_sorted.index.values
    sigPValInds = indexVals[:5]
    return tpValDf, tpValDf_sorted, sigPValInds
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
def showMePCAFig(DataSlice, targetCol):
    ''' This function displays a figure derived from PCA '''

    from sklearn.preprocessing import StandardScaler

    # Standardizing the features
    x = StandardScaler().fit_transform(DataSlice)

    from sklearn.decomposition import PCA

    pca = PCA(n_components=2)
    principalComponents = pca.fit_transform(x)
    principalDf = pd.DataFrame(data=principalComponents, columns=['PC1', 'PC2'])

    targetCol_df = pd.DataFrame(data=targetCol, columns=['Diagnosis'])

    finalDf = pd.concat([principalDf, targetCol_df], axis = 1)

    # Plotting the PCA scatterplot using Seaborn
    sns.set()
    ax = sns.relplot(x='PC1',y='PC2',data=finalDf,hue='Diagnosis',palette='Set1')
    plt.show()
    return
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
def showMeViolinPlts(targetCol, sigPValInds, DataSlice, boolean, number):
    ''' This function creates violin plots using 4 main variables:
    targetCol, sigPValInds, DataSlice, boolean
    # boolean == 0: when looking at FeatSlices
    # boolean == 1: when looking at ROISlices '''

    # Create an index for the subplot
    n = 1;
    fig = plt.figure()

    # Copy pasted (code enables custom selection of rows)
    a = np.where(targetCol == 0)[0]
    b = np.where(targetCol == 1)[0]

    # Z-score the data
    DataSlice = pd.DataFrame(data=DataSlice)
    DataSlice = DataSlice.apply(zscore)

    for i in sigPValInds:
        # Obtain control feature i and SCZ feature i (all the rows)
        # from the ith column of the z-scored feature matrix
        cf_i = DataSlice.iloc[a,i-1]
        sf_i = DataSlice.iloc[b,i-1]

        # Convert into a dataframe
        df_feat_i = pd.DataFrame({'SCZ':sf_i,'Control':cf_i})

        # Violin plots
        ax = fig.add_subplot(2,3,n)
        ax = sns.violinplot(data=df_feat_i, order=['SCZ','Control'])
        plt.xlabel('Diagnosis')

        if boolean == 0:
            ylabel = 'ROI ' + str(i)
        elif boolean == 1:
            ylabel = 'Feature ' + str(i)
        plt.ylabel(ylabel)

        # Increment index
        n += 1;

    if boolean == 0:
         plt.suptitle('Top 5 ROI in Feature ' + str(number), fontsize=12.5, y=1.05)
    elif boolean == 1:
        plt.suptitle('Top 5 Features in ROI ' + str(number), fontsize=12.5, y=1.05)
    plt.tight_layout()
    plt.show()
    return
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
def showMeROIAccPlot(ROIs, tsData, indices2KeepMat, targetCol, dispFigs):
    ''' This function displays the region accuracy plot and the top 5 regions
    having the highest accuracies
    If 0 is given as an input, this function suppresses any printed messages
    and plots but returns the mean regional accuracy and error '''

    # Initialise a few variables
    regions = np.zeros([ROIs])
    region_acc = np.zeros([ROIs])
    roiErr = np.zeros([ROIs])

    for n in range(1, ROIs+1):

        # For each of the 1 to n regions, show me how good all of the features are at predicting
        # whether the subject has SCZ or not
        tsDataSlice = tsData.loc[n,indices2KeepMat,:]
        tsDataSlice_zscored = tsDataSlice.apply(zscore)

        ## Print the top 5 features for each region being analysed
        # top5Feats = getTPVals(targetCol, tsDataSlice)[2]
        # print(top5Feats)

        # Assign the data to variables
        X = tsDataSlice_zscored
        y = np.ravel(targetCol)

        # Denotes how good all the features in a given region are at predicting if a subject has SCZ
        avgScore = get10FoldCVScore(X,y).mean()
        avgSTD = get10FoldCVScore(X,y).std()
        regions[n-1] = n
        region_acc[n-1] = avgScore
        roiErr[n-1] = avgSTD

    df = pd.DataFrame({'Region':regions,'% Accuracy':region_acc,'ROI Error':roiErr})
    df['Region'] = df.Region.astype(int)

    df_sorted = df.sort_values(by='% Accuracy',ascending=False)
    df_sorted = df_sorted.set_index('Region')

    if dispFigs == 1:
        print('')
        print(df_sorted[:5].to_string())
        print('')

        print('Mean Accuracy (across all regions) = ' +
        "{0:.2f}".format(df['% Accuracy'].mean()) + '%')
        print('')

        print('Mean Error (across all regions) = ' +
        "{0:.2f}".format(df['ROI Error'].mean()) + '%')
        print('')

        sns.distplot(region_acc, bins='auto', kde=False)
        plt.xlabel('Classification Accuracy (%)')
        plt.ylabel('No. of Regions')
        plt.show()
        return

    elif dispFigs == 0:
        meanROIAcc = '{0:.2f}'.format(df['% Accuracy'].mean())
        meanROIError = '{0:.2f}'.format(df['ROI Error'].mean())
        return meanROIAcc, meanROIError
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
def showMeFeatAccPlot(element, subPath, PyFeatList, indices2KeepMat, targetCol, dispFigs):
    ''' This function displays the feature accuracy plot and the top 5 features
    having the highest accuracies
    If 0 is given as an input, this function suppresses any printed messages
    and plots but returns the mean feature accuracy and error '''

    tsData, ROIs, subjects, feats, featList = addIndices(element, subPath, PyFeatList)

    # Initialise a few variables
    featNo = np.zeros([feats])
    feat_acc = np.zeros([feats])
    featErr = np.zeros([feats])

    for n in range(1, feats+1):

        # For each of the 1 to 22 features, show me how good all of the ROIs are at predicting
        # whether the subject has SCZ or not
        tsDataSlice = getFeatSlice(ROIs, subjects, tsData, featList[n-1], indices2KeepMat)
        tsDataSlice_zscored = tsDataSlice.apply(zscore)

        # Assign the data to variables
        X = tsDataSlice_zscored
        y = np.ravel(targetCol)

        avgScore = get10FoldCVScore(X,y).mean()
        avgSTD = get10FoldCVScore(X,y).std()
        featNo[n-1] = n
        feat_acc[n-1] = avgScore
        featErr[n-1] = avgSTD

    df = pd.DataFrame({'Feature':featNo,'% Accuracy':feat_acc,'Feat Error':featErr})
    df['Feature'] = df.Feature.astype(int)

    df_sorted = df.sort_values(by='% Accuracy',ascending=False)
    df_sorted = df_sorted.set_index('Feature')

    if dispFigs == 1:
        print(df_sorted[:5].to_string())
        print('')

        print('Mean Accuracy (across all features) = ' +
        "{0:.2f}".format(df['% Accuracy'].mean()) + '%')
        print('')

        print('Mean Error (across all features) = ' +
        "{0:.2f}".format(df['Feat Error'].mean()) + '%')
        print('')

        sns.distplot(feat_acc, bins='auto', kde=False)
        plt.xlabel('Classification Accuracy (%)')
        plt.ylabel('No. of Features')
        plt.show()
        return

    elif dispFigs == 0:
        meanFeatAcc = '{0:.2f}'.format(df['% Accuracy'].mean())
        meanFeatError = '{0:.2f}'.format(df['Feat Error'].mean())
        return meanFeatAcc, meanFeatError
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
def Reg_by_Reg_Anal(ROI, tsData, targetCol, ROIs, indices2KeepMat, regAccOnly, dispFigs):
    ''' This function takes in several inputs that have been calculated and
    outputs several graphs pertaining to the analysis of the selected region -
    It provides a Region-by-Region analysis, and the graphs displayed can be
    suppressed by setting the variable 'dispFigs' as false '''

    # Acquire ROI slice
    ROISlice = tsData.loc[ROI,indices2KeepMat,:]

    # Assign the data to variables
    DataSlice = ROISlice
    DataSlice_zscored = DataSlice.apply(zscore)

    X = DataSlice_zscored
    y = np.ravel(targetCol)

    # Perform 10-fold CV
    scores = get10FoldCVScore(X,y)
    avgScore = '{0:.2f}'.format(scores.mean())
    avgSTD = '{0:.2f}'.format(scores.std())

    if regAccOnly == True:
        return avgScore, avgSTD

    if dispFigs == True:

        # Print scores
        print('Analysis of Region ' + str(ROI) + ':')
        print('')

        print('10-fold CV scores as a percentage: ' + str(scores))
        print('')

        # Mean 10-fold CV score with an error of 1 std dev
        print("Accuracy as a percentage: " + avgScore + " +/- " + avgSTD)

        # Store the first five indices of the features with the most significant p-values (the third output)
        tpValDf, tpValDf_sorted, sigPValInds = getTPVals(targetCol, DataSlice)

        # Show me the PCA figure
        showMePCAFig(DataSlice, targetCol)

        # Show me the top five features as violin plots in the ROI being analysed
        showMeViolinPlts(targetCol, sigPValInds, DataSlice, 1, ROI)

        # Show me the average classification accuracy of all the features in a
        # given region, and then take the average of all the regions
        showMeROIAccPlot(ROIs, tsData, indices2KeepMat, targetCol, dispFigs)
        return

    elif dispFigs == False:
        meanROIAcc, meanROIError = showMeROIAccPlot(ROIs, tsData, indices2KeepMat, targetCol, dispFigs)
        return avgScore, avgSTD, meanROIAcc, meanROIError
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
def Feat_by_Feat_Anal(feature, featureName, element, subPath, PyFeatList,
indices2KeepMat, targetCol, featAccOnly, dispFigs):

    tsData, ROIs, subjects, feats, featList = addIndices(element, subPath, PyFeatList)

    # Acquire feature slice
    featSlice = getFeatSlice(ROIs, subjects, tsData, featureName, indices2KeepMat)

    # Assign the data to variables
    DataSlice = featSlice
    DataSlice_zscored = DataSlice.apply(zscore)

    X = DataSlice_zscored
    y = np.ravel(targetCol)

    # Perform 10-fold CV
    scores = get10FoldCVScore(X,y)
    avgScore = '{0:.2f}'.format(scores.mean())
    avgSTD = '{0:.2f}'.format(scores.std())

    if featAccOnly == True:
        return avgScore, avgSTD

    if dispFigs == True:

        # Print scores
        print('Analysis of Feature ' + str(feature) + ':')
        print('')

        print('10-fold CV scores as a percentage: ' + str(scores))
        print('')

        # Mean 10-fold CV score with an error of 1 std dev
        print("Accuracy as a percentage: " + avgScore + " +/- " + avgSTD)

        # Store the first five indices of the ROIs with the most significant p-values (the third output)
        tpValDf, tpValDf_sorted, sigPValInds = getTPVals(targetCol, DataSlice)

        # Show me the PCA figure
        showMePCAFig(DataSlice, targetCol)

        # Show me the top five ROIs as violin plots in the feature being analysed
        showMeViolinPlts(targetCol, sigPValInds, DataSlice, 0, feature)

        # Show me the average classification accuracy of all the ROIs in a
        # given feature, and then take the average of all the features
        showMeFeatAccPlot(element, subPath, PyFeatList, indices2KeepMat, targetCol, dispFigs)
        return

    elif dispFigs == False:
        meanFeatAcc, meanFeatError = showMeFeatAccPlot(element, subPath, PyFeatList,
        indices2KeepMat, targetCol, dispFigs)
        return avgScore, avgSTD, meanFeatAcc, meanFeatError
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
def giveMeFDvBalancedAcc(fdArray):

    fig, ax1 = plt.subplots()

    AvgROIError = fdArray.iloc[:,6]
    AvgFeatError = fdArray.iloc[:,8]

    fd = fdArray.iloc[:,0]
    SCZ2Control = fdArray.iloc[:,4]
    AvgROIAcc = fdArray.iloc[:,5]
    AvgFeatAcc = fdArray.iloc[:,7]

    ax1.plot(fd, AvgFeatAcc,'b-')
    # ax1.errorbar(fd, AvgFeatAcc, yerr=AvgFeatError, fmt='-o', capsize=5)

    ax1.plot(fd, AvgROIAcc, 'r-')
    # ax1.errorbar(fd, AvgROIAcc, yerr=AvgROIError, fmt='-o', capsize=5)

    plt.xlim(max(fd)+0.01, min(fd)-0.01)
    ax1.set_xlabel('fd (cm)')
    ax1.set_ylabel('Balanced Acc (%)')

    ax2 = ax1.twinx()
    ax2.plot(fd, SCZ2Control,'g-')
    ax2.set_ylabel('SCZ:Control', color='g')
    ax2.tick_params('y', colors='g')
    plt.grid(b=None)

    ax1.legend(['Avg Feat Acc','Avg ROI Acc'],loc=3)
    ax2.legend(['SCZ:Control'],loc=2)
    fig.tight_layout()
    plt.show()
    return
#-------------------------------------------------------------------------------
