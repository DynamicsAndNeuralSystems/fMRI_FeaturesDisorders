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
def addIndices(element_txt,subPath,PyFeatList_txt):
    ''' This function adds a multi-level index to the feature matrix data that
        is read in as a text file, and returns the matrix
        The function takes two other inputs: a .txt file of the feature names
        as well as the subject file path for the given dataset '''

    from numpy import genfromtxt
    tsData = genfromtxt(element_txt, delimiter=',')
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

    featList = [lines.rstrip('\n') for lines in open(PyFeatList_txt)]
    tsData = pd.DataFrame(data=tsData, index=index, columns=featList)

    return tsData, ROIs, subjects, feats
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
def getTargetCol(TS_path_names):
    ''' This function, when called, prompts the user to input a string
    which can be used to uniquely identify any SCZ data file
    It returns a binary target column (0 = control, 1 = SCZ) '''

    # string = input('Enter a string found in all filenames that only contain SCZ data: ')
    # In this case '-5' should do
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
def getROISlice(Orig_TS_path_names, tsData, ROI, indices2Keep):
    ''' This function takes a ROI as an input and gets the relevant slice
    (or rows) from tsData '''

#     ROI = input('Which ROI is being analysed?: ')

    # Initialise a boolean
    validROI = False

    while validROI == False:

        # Rows per region of interest = len(TS_path_names)
        rowsPerROI = len(Orig_TS_path_names)
        [rows,cols] = tsData.shape

        # Get the rows corresponding to the mth region of interest
        ROISlice = tsData.iloc[((int(ROI)-1)*rowsPerROI):(int(ROI)*rowsPerROI),:]
        ROISlice = ROISlice.iloc[indices2Keep,:]

        # Assigning min and max values for input 'm'
        minROI = 1
        maxROI = int(rows/rowsPerROI)

        if minROI <= int(ROI) <= maxROI:
            validROI == True
            return ROISlice, ROI, maxROI
#         else:
#             print('Error: Please select a ROI between ' + str(minROI) + ' and ' + str(maxROI))
#             print('')
#             ROI = input('Which ROI is being analysed?: ')
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
def getFeatSlice(featMat3D, feature):
    ''' This function returns a slice of all the data related to a particular
    feature, which it receives as an integer input '''

    if 0 < feature <= 22:
        FeatSlice = featMat3D[:,feature-1,:]
        return FeatSlice
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
def get10FoldCVScore(X,y):
    ''' This function returns a 10-fold CV score after balancing the classes '''

    # Import the support vector classifier and balance the classes
    from sklearn.svm import SVC
    svclassifier = SVC(kernel='linear')

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

        scores[i] = '{0:.1f}'.format(balanced_accuracy_score(y_test, y_pred)*100)
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

    # Sort the data (including the indices) in descending order by MAGNITUDE
    tpValDf_sorted = tpValDf.abs().sort_values(by='t-value',ascending=False)

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

    # boolean = int(input('Is DataSlice = FeatSlice or ROISlice? Enter 0 or 1 respectively: '))
    # number = input('Which Feature / ROI is being analysed?: ')

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
        cf_i = DataSlice.iloc[a,i]
        sf_i = DataSlice.iloc[b,i]

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
def showMeROIAccPlot(maxROI, Orig_TS_path_names, tsData, targetCol, boolean, indices2Keep):
    ''' This function displays the region accuracy plot and the top 5 regions
    having the highest accuracies
    If 0 is given as an input, this function suppresses any printed messages
    and plots but returns the mean regional accuracy '''

    # Initialise a few variables
    regions = np.zeros([maxROI])
    region_acc = np.zeros([maxROI])
    roiErr = np.zeros([maxROI])

    for n in range(1, maxROI+1):

        tsDataSlice = getROISlice(Orig_TS_path_names, tsData, n, indices2Keep)[0]
        tsDataSlice_zscored = tsDataSlice.apply(zscore)

        # Assign the data to variables
        X = tsDataSlice_zscored
        y = np.ravel(targetCol)

        avgScore = get10FoldCVScore(X,y).mean()
        avgSTD = get10FoldCVScore(X,y).std()
        regions[n-1] = n
        region_acc[n-1] = avgScore
        roiErr[n-1] = avgSTD

    df = pd.DataFrame({'Region':regions,'% Accuracy':region_acc,'ROI Error':roiErr})
    df['Region'] = df.Region.astype(int)

    df_sorted = df.sort_values(by='% Accuracy',ascending=False)

    if boolean == 1:
        print(df_sorted[:5].to_string(index=False))
        print('')

        print('Mean Accuracy (across all regions) = ' +
        "{0:.2f}".format(df['% Accuracy'].mean()) + '%')
        print('')

        print('Mean Error (across all regions) = ' +
        "{0:.2f}".format(df['ROI Error'].mean()) + '%')
        print('')

        plt.hist(region_acc, bins='auto')
        plt.xlabel('Classification Accuracy (%)')
        plt.ylabel('No. of Regions')
        plt.show()
        return
    elif boolean == 0:
        meanROIAcc = '{0:.2f}'.format(df['% Accuracy'].mean())
        AvgROIError = '{0:.2f}'.format(df['ROI Error'].mean())
        return meanROIAcc, AvgROIError
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
def showMeFeatAccPlot(featMat3D, targetCol, boolean):
    ''' This function displays the feature accuracy plot and the top 5 features
    having the highest accuracies
    If 0 is given as an input, this function suppresses any printed messages
    and plots but returns the mean feature accuracy '''

    # Initialise a few variables
    feats = np.zeros(22)
    feat_acc = np.zeros(22)
    featErr = np.zeros([22])

    for n in range(1, 23):

        tsDataSlice = pd.DataFrame(getFeatSlice(featMat3D, n))
        tsDataSlice_zscored = tsDataSlice.apply(zscore)

        # Assign the data to variables
        X = tsDataSlice_zscored
        y = np.ravel(targetCol)

        avgScore = get10FoldCVScore(X,y).mean()
        avgSTD = get10FoldCVScore(X,y).std()
        feats[n-1] = n
        feat_acc[n-1] = avgScore
        featErr[n-1] = avgSTD

    df = pd.DataFrame({'Feature':feats,'% Accuracy':feat_acc,'Feat Error':featErr})
    df['Feature'] = df.Feature.astype(int)

    df_sorted = df.sort_values(by='% Accuracy',ascending=False)

    if boolean == 1:
        print(df_sorted[:5].to_string(index=False))
        print('')

        print('Mean Accuracy (across all features) = ' +
        "{0:.2f}".format(df['% Accuracy'].mean()) + '%')
        print('')

        print('Mean Error (across all features) = ' +
        "{0:.2f}".format(df['Feat Error'].mean()) + '%')
        print('')

        plt.hist(feat_acc, bins='auto')
        plt.xlabel('Classification Accuracy (%)')
        plt.ylabel('No. of Features')
        plt.show()
        return
    elif boolean == 0:
        meanFeatAcc = '{0:.2f}'.format(df['% Accuracy'].mean())
        AvgFeatError = '{0:.2f}'.format(df['Feat Error'].mean())
        return meanFeatAcc, AvgFeatError, df_sorted
#-------------------------------------------------------------------------------
