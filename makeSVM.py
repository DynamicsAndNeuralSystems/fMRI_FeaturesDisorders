#!/usr/local/bin/python3
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('classic')
import seaborn as sns
import sklearn
from scipy.stats import stats,ttest_ind

# Import and store the feature matrix in the variable tsData
path = 'featureMatrixPy.txt'
tsData = pd.read_csv(path,header=None);
# print(tsData)

"NEED TO GENERALISE FOR ANY CASE"
[rows, cols] = tsData.shape

# bool = input('Is the first half of the data set the control set? y/n ')

# Create a 'target' column where rows 0-99 have the value 1 (indicating a seizure)
# and rows 100-199 have the value 0 (indicating no seizure)
zeros = np.zeros(100, dtype=int)
ones = np.ones(100, dtype=int)
targetCol = np.hstack((ones, zeros))
targetCol = np.reshape(targetCol,(200,1))

# Need to z-score the feature matrix and save it back into the variable, tsData
# We are looking down each column and calculating the z-score for each feature
from scipy.stats import zscore
tsData_zscored = tsData.apply(zscore)
# print(tsData)

# Assign data to variables
X = tsData_zscored
y = np.ravel(targetCol)

# Split the data into training and test sets
from sklearn.model_selection import train_test_split
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.10)

# Import the support vector classifier and find the best fit using the training data
from sklearn.svm import SVC
svclassifier = SVC(kernel='linear')
svclassifier.fit(X_train, y_train)

# Using the SVC derived from fitting to the training set, classify the test set
# and store this in a variable
y_pred = svclassifier.predict(X_test)

# Compare the predictions made by the SVC to the correct classifications
# Display the results using a confusion matrix and the classification accuracy
from sklearn.metrics import accuracy_score, confusion_matrix
print(confusion_matrix(y_test, y_pred))
print('')
print('The classification accuracy of the SVM for a single random split is ' +
str(accuracy_score(y_test, y_pred) * 100) + '%')

# Calculate specificity (true negatives) and sensitivity (true positives)
confMat = confusion_matrix(y_test, y_pred)
specificity = (confMat[0,0] / (confMat[0,0] + confMat[1,0])) * 100
sensitivity = (confMat[1,1] / (confMat[1,1] + confMat[0,1])) * 100
print('The specificity is ' + str("%.1f" % specificity) + '% and the sensitivity is '
+ str("%.1f" % sensitivity) + '%')
print('')

"NEED TO CALCULATE THE AVERAGE SPECIFICITY AND SENSITIVITY"
# Performing a 10-fold validation using cross_val_score
from sklearn.model_selection import cross_val_score
scores = cross_val_score(svclassifier, X, y, cv=10) * 100
scores = scores.astype(int)

# Print scores
print('10-fold CV scores as a percentage: ' + str(scores))
print('')

avgCVScore = np.mean(scores)
print('The average 10-fold CV score is ' + str("%.1f" % avgCVScore) + '%')
print('')

# Compute the t-value (from a two-tail t-test) and the p-value
# Store these two values (t-value, then p-value) in each row, 22 in total for each feature

# Initialise the array and assign its shape
tpValArray = np.zeros([22, 2])
[rows, cols] = tpValArray.shape

"NEED TO GENERALISE FOR ANY CASE"
# Loop through the array and store the t and p values
for i in range(rows):

    # Calculate the t and p values by inputting the two halves of each of the 22 columns
    # of the normalised data into the ttest functions
    # Store the statistics in the variable, tpVal (which changes on each iteration of the outer loop)
    setE_FeatureCol = tsData.iloc[0:100,i]
    setA_FeatureCol = tsData.iloc[100:200,i]
    tpVal = stats.ttest_ind(setA_FeatureCol, setE_FeatureCol)

    for j in range(cols):

        # Store the values into each column
        tpValArray[i,j] = tpVal[j]

# Since it is a two-tailed t-test, need to multiply the p-values by two (second column)
tpValArray[:,1] = tpValArray[:,1] * 2

# Formatting the tpValArray
tpValDf = pd.DataFrame(data=tpValArray, columns=['t-value', 'p-value'])
tpValDf.index.name = 'Feature i'
print(tpValDf)
print('')

# Sort the data (including the indices) in descending order by MAGNITUDE
tpValDf_sorted = tpValDf.abs().sort_values(by='t-value',ascending=False)
print(tpValDf_sorted)
print('')

# Store the first five indices of the sorted dataframe - will need to use these
# indices to access the relevant feature columns in X, the feature matrix
indexVals = tpValDf_sorted.index.values
signifTVals = indexVals[:5]
print(signifTVals)
print('')

#-------------------------------------------------------------------------------
# PCA
from sklearn.preprocessing import StandardScaler

# Standardizing the features
x = StandardScaler().fit_transform(tsData)

from sklearn.decomposition import PCA
pca = PCA(n_components=2)
principalComponents = pca.fit_transform(x)
principalDf = pd.DataFrame(data=principalComponents
             , columns=['PC1', 'PC2'])

targetCol_df = pd.DataFrame(data=targetCol, columns=['target'])

finalDf = pd.concat([principalDf, targetCol_df], axis = 1)

print(finalDf)

# Plotting the PCA scatterplot using Seaborn
sns.set()
ax = sns.relplot(x='PC1', y='PC2', data=finalDf, hue='target',palette='Set1')
plt.show()
#-------------------------------------------------------------------------------
# Create an index for the subplot
n = 1;

fig = plt.figure()

for i in signifTVals:
    # Obtaining Seizure_Feature_i & Healthy_Feature_i
    # (all the rows) from the ith column of the feature matrix
    sf_i = X.iloc[0:100,i]
    hf_i = X.iloc[100:200,i]

    # Stacking columns side by side
    feat_i = np.column_stack((sf_i,hf_i))

    # Making the numpy array into a dataframe
    df_feat_i = pd.DataFrame(data=feat_i, columns=['Epileptic', 'Healthy'])
    # print(df_feat_i)

    # Violin plots
    ax = fig.add_subplot(2,3,n)
    ax = sns.violinplot(data=df_feat_i, order=["Epileptic", "Healthy"])
    plt.xlabel('Diagnosis')
    ylabel = 'Feature ' + str(i)
    plt.ylabel(ylabel)

    # Increment index
    n += 1;

plt.tight_layout()
plt.show()

# pVal = tpValDf.iloc[:,1]
#
# from scipy.interpolate import spline
# x = np.arange(22)
# print(x)
# x_smooth = np.linspace(x.min(),x.max(),300)
# pVal_smooth = spline(x,pVal,x_smooth)
#
# plt.plot(x_smooth,pVal_smooth)
# plt.show()
