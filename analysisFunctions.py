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
def removePathNames(threshold_fd, TS_path_names):
    ''' This function modifies the entries within TS_path_names based on the
    threshold fd value '''

    filePath = '/Users/AV/Dropbox/UCLA/movementData/fdAvgs.txt'
    fdAvgs = pd.read_csv(filePath,header=None);

    indices = np.where(fdAvgs > threshold_fd)[0]

    for i in sorted(indices, reverse=True):
        del TS_path_names[i]
    return TS_path_names
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
def getROISlice(TS_path_names, tsData, ROI):
    ''' This function takes a ROI as an input and gets the relevant slice
    (or rows) from tsData '''

#     ROI = input('Which ROI is being analysed?: ')

    # Initialise a boolean
    validROI = False

    while validROI == False:

        # Rows per region of interest = len(TS_path_names)
        rowsPerROI = len(TS_path_names)
        [rows,cols] = tsData.shape

        # Get the rows corresponding to the mth region of interest
        Slice = tsData.iloc[((int(ROI)-1)*rowsPerROI):(int(ROI)*rowsPerROI),:]

        # Assigning min and max values for input 'm'
        minROI = 1
        maxROI = int(rows/rowsPerROI)

        if minROI <= int(ROI) <= maxROI:
            validROI == True
            return Slice, ROI, maxROI
#         else:
#             print('Error: Please select a ROI between ' + str(minROI) + ' and ' + str(maxROI))
#             print('')
#             ROI = input('Which ROI is being analysed?: ')
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
def getFeatSlice(featMat3D, feat):
    ''' This function returns a slice of all the data related to a particular
    feature, which it receives as an integer input '''

    Slice = featMat3D[:,abs(feat-1),:]
    return Slice
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
def get10FoldCVScore(X,y):
    ''' This function returns a 10-fold CV score after balancing the classes '''

    # Import the support vector classifier and balance the classes
    from sklearn.svm import SVC
    svclassifier = SVC(kernel='linear', class_weight = 'balanced')

    # Split the data into training and test sets
    from sklearn.model_selection import KFold
    kf = KFold(n_splits=10)

    # Import accuracy score
    from sklearn.metrics import accuracy_score

    # Initialise a few variables
    scores = np.zeros(10)
    i = 0

    for train_index, test_index in kf.split(X):

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

        scores[i] = '{0:.1f}'.format(accuracy_score(y_test, y_pred)*100)
#         print('Acc % = ', scores[i])
#         print('')

        # Increment index
        i += 1
    return scores
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
def getTPVals(targetCol, tsDataSlice):
    ''' This function computes the t-values (from a two-tail t-test) and
    the p-values
    It stores these two values (t-value, then p-value) in each row,
    22 in total for each feature
    It returns a few outputs, signifTVals being the main one '''

    from scipy.stats import stats,ttest_ind

    # Initialise the array and assign its shape
    tpValArray = np.zeros([22, 2])
    [rows, cols] = tpValArray.shape

    # Make a custom range of rows which can be used to distinguish the control data from the SCZ data
    # Find all instances of 0 in the target column (which would indicate control data)
    # and store the indices in an array
    a = np.where(targetCol == 0)[0]
    aStart = a[0]
    aEnd = a[-1] + 1

    # Find all instances of 1 in the target column (which would indicate SCZ data)
    # and store the indices in an array
    b = np.where(targetCol == 1)[0]
    bStart = b[0]
    bEnd = b[-1] + 1

    # Loop through the array and store the t and p values
    for i in range(rows):

        # Calculate the t and p values by inputting the two halves of each of the 22 columns
        # of the normalised data into the ttest functions
        # Store the statistics in the variable, tpVal (which changes on each iteration of the outer loop)
        controlFeatCol = tsDataSlice.iloc[aStart:aEnd,i]
        SCZFeatCol = tsDataSlice.iloc[bStart:bEnd,i]
        tpVal = stats.ttest_ind(controlFeatCol, SCZFeatCol)

        for j in range(cols):

            # Store the values into each column
            tpValArray[i,j] = tpVal[j]

    # Since it is a two-tailed t-test, need to multiply the p-values by two (second column)
    tpValArray[:,1] = tpValArray[:,1] * 2

    # Formatting the tpValArray
    tpValDf = pd.DataFrame(data=tpValArray, columns=['t-value', 'p-value'])
    tpValDf.index.name = 'Feature i'

    # Sort the data (including the indices) in descending order by MAGNITUDE
    tpValDf_sorted = tpValDf.abs().sort_values(by='t-value',ascending=False)

    # Store the first five indices of the sorted dataframe - will need to use these
    # indices to access the relevant feature columns in X, the feature matrix
    indexVals = tpValDf_sorted.index.values
    signifTVals = indexVals[:5]

    return tpValDf, tpValDf_sorted, signifTVals
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
def showMePCAFig(tsDataSlice, targetCol):
    ''' This function displays a figure derived from PCA '''

    from sklearn.preprocessing import StandardScaler

    # Standardizing the features
    x = StandardScaler().fit_transform(tsDataSlice)

    from sklearn.decomposition import PCA
    pca = PCA(n_components=2)
    principalComponents = pca.fit_transform(x)
    principalDf = pd.DataFrame(data=principalComponents
                 , columns=['PC1', 'PC2'])

    targetCol_df = pd.DataFrame(data=targetCol, columns=['target'])

    finalDf = pd.concat([principalDf, targetCol_df], axis = 1)

    # Plotting the PCA scatterplot using Seaborn
    sns.set()
    ax = sns.relplot(x='PC1',y='PC2',data=finalDf,hue='target',palette='Set1')
    plt.show()
    return
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
def showMeViolinPlts(targetCol, signifTVals, X, ROI):
    ''' This function creates violin plots using 4 main variables:
    targetCol, signifTVals, X, ROI - probably can be modified a little to
    plot the features as well as the regions '''

    # Create an index for the subplot
    n = 1;

    fig = plt.figure()

    # Copy pasted (code enables custom selection of rows)
    a = np.where(targetCol == 0)[0]
    aStart = a[0]
    aEnd = a[-1] + 1

    b = np.where(targetCol == 1)[0]
    bStart = b[0]
    bEnd = b[-1] + 1

    for i in signifTVals:
        # Obtain control feature i and SCZ feature i (all the rows)
        # from the ith column of the z-scored feature matrix
        cf_i = X.iloc[aStart:aEnd,i]
        sf_i = X.iloc[bStart:bEnd,i]

        # Convert into a dataframe
        df_feat_i = pd.DataFrame({'SCZ':sf_i,'Control':cf_i})

        # Violin plots
        ax = fig.add_subplot(2,3,n)
        ax = sns.violinplot(data=df_feat_i, order=['SCZ','Control'])
        plt.xlabel('Diagnosis')
        ylabel = 'Feature ' + str(i)
        plt.ylabel(ylabel)

        # Increment index
        n += 1;

    plt.suptitle('Top 5 Features in ROI ' + str(ROI), fontsize=12.5, y=1.05)
    plt.tight_layout()
    plt.show()
    return
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
def showMeRegAccPlot(maxROI, TS_path_names, tsData, targetCol, bool):
    ''' This function displays the region accuracy plot and the top 5 regions
    having the highest accuracies
    If 0 is given as an input, this function suppresses any printed messages
    and plots but returns the mean regional accuracy '''

    # Initialise a few variables
    regions = np.zeros([maxROI])
    region_acc = np.zeros([maxROI])

    for n in range(1, maxROI+1):

        tsDataSlice = getROISlice(TS_path_names, tsData, n)[0]
        tsDataSlice_zscored = tsDataSlice.apply(zscore)

        # Assign the data to variables
        X = tsDataSlice_zscored
        y = np.ravel(targetCol)

        avgScore = get10FoldCVScore(X,y).mean()
        regions[n-1] = n
        region_acc[n-1] = avgScore

    df = pd.DataFrame({'Region':regions,'% Accuracy':region_acc})
    df['Region'] = df.Region.astype(int)

    df_sorted = df.sort_values(by='% Accuracy',ascending=False)

    if bool == 1:
        print(df_sorted[:5].to_string(index=False))
        print('')

        print('Mean Accuracy (across all regions) = ' +
        "{0:.2f}".format(df['% Accuracy'].mean()) + '%')
        print('')

        plt.hist(region_acc, bins='auto')
        plt.xlabel('Classification Accuracy (%)')
        plt.ylabel('No. of Regions')
        plt.show()
        return
    elif bool == 0:
        meanRegAcc = '{0:.2f}'.format(df['% Accuracy'].mean())
        return meanRegAcc
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
def showMeFeatAccPlot(featMat3D, targetCol, bool):
    ''' This function displays the feature accuracy plot and the top 5 features
    having the highest accuracies
    If 0 is given as an input, this function suppresses any printed messages
    and plots but returns the mean feature accuracy '''

    # Initialise a few variables
    feats = np.zeros(22)
    feat_acc = np.zeros(22)

    for n in range(1, 23):

        tsDataSlice = pd.DataFrame(getFeatSlice(featMat3D, n))
        tsDataSlice_zscored = tsDataSlice.apply(zscore)

        # Assign the data to variables
        X = tsDataSlice_zscored
        y = np.ravel(targetCol)

        avgScore = get10FoldCVScore(X,y).mean()
        feats[n-1] = n
        feat_acc[n-1] = avgScore

    df = pd.DataFrame({'Feature':feats,'% Accuracy':feat_acc})
    df['Feature'] = df.Feature.astype(int)

    df_sorted = df.sort_values(by='% Accuracy',ascending=False)

    if bool == 1:
        print(df_sorted[:5].to_string(index=False))
        print('')

        print('Mean Accuracy (across all features) = ' +
        "{0:.2f}".format(df['% Accuracy'].mean()) + '%')
        print('')

        plt.hist(feat_acc, bins='auto')
        plt.xlabel('Classification Accuracy (%)')
        plt.ylabel('No. of Features')
        plt.show()
        return
    elif bool == 0:
        meanFeatAcc = '{0:.2f}'.format(df['% Accuracy'].mean())
        return meanFeatAcc
#-------------------------------------------------------------------------------
