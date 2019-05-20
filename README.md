# _catch22_ for cfg.mat files
There are three main files, which will:
* Firstly, read in a _cfg.mat_ file and output an _'elementnDATASET.txt'_ file, which is the format that the time series data will be stored in
* Secondly, clearly define all of the analysis funtions that will be used
* And finally, call upon the analysis functions and save the outputs

The main functions in the analysis (along with the respective lines of code) are as follows:

#### Generating the Time Series Data text file
First, define a path and sort the subject filenames:
```python
path = '/Users/AV/Dropbox/COBRE/cfgData/'
TS_path_names = sorted(glob.glob(path + '*.mat'))
```

Next, extract the time series data from the cfg.mat file and store it in a variable:
```python
gettsData(path,element)
```

For a given time series, calculate the 1 x 22 feature vector using _catch22_:
```python
giveMeFeatureVector(tsDataList)
```

For a given region of interest, ROI, (a column of the tsData matrix), calculate the feature vector for _all_ subjects, and store this matrix in a variable called 'smallFeatMat':
```python
giveMeSmallFeatMat(element,roi)
```

Continue getting the smallFeatMats for each ROI, stack all of these matrices and store it in a variable called 'bigFeatMat':
```python
giveMeBigFeatMat(element)
```

Note, the function _giveMeBigFeatMat_ utilises all the above functions.
Call the _giveMeBigFeatMat_ function and output the tsData text file:
```python
# Which element of 'roiTS' should be analysed?
x, element = giveMeBigFeatMat(1)

# Save file
np.savetxt('element1_COBRE' + str(element) + '.txt',x,delimiter=',')
```

#### The Analysis Functions

Since the effect of fd on the classification accuracy needed to be observed, the first few lines of code (as well as the first defined function) will address this aim.
An fd array has been created (using the fdAvgs) which will be used within a for loop:
```python
fdArray = np.linspace(minFd, maxFd, n)
for threshold_fd in fdArray:
```

Next, the subject files need to be filtered based on the subjects' fdAvgs:
```python
removePathNames(filePath, threshold_fd, TS_path_names)

# Store the new TS_path_names and the indices of the subject files that are being kept
TS_path_names, indices2Keep = af.removePathNames(filePath, threshold_fd, TS_path_names)
```

Create a target column **(haven't made this work for the general case yet!)**:
```python
targetCol = af.getTargetCol(TS_path_names)
```

Define a function that can take 'ROI slices' from the bigFeatMat.
Recall that bigFeatMat is a stack of all the smallFeatMats, with each smallFeatMat containing the feature vectors for all subjects for a particular ROI:
```python
getROISlice(Orig_TS_path_names, tsData, ROI, indices2Keep)
```

Reshape the feature data into a 3D matrix with dimensions defined by the number of subjects x the number of features x the number of ROIs:
```python
featMat3D = np.zeros((len(TS_path_names),22,maxROI))

# Loop through each slice and store it in the matrix
for i in range(1, maxROI+1):
    featMat3D[:,:,i-1] = af.getROISlice(Orig_TS_path_names, tsData, i, indices2Keep)[0]
```

Define another function that can access the 'Feature slices' from the 3D feature matrix, _featMat3D_:
```python
getFeatSlice(featMat3D, feature)
```

Calculate the 10-Fold CV using StratifiedKFold and balanced_accuracy_score:
```python
# This function needs to be given either a 'ROISlice' or a 'FeatSlice' (X) and the target column (y)
# Return an array of balanced classification accuracies, one for each of the 10 folds
get10FoldCVScore(X,y)
```

Calculate the t- and p-test values and return the values that are significant:
```python
# This function computes and returns the t- and p-test values for each of the features or ROIs
getTPVals(targetCol, DataSlice)
```

Generate a bunch of different plots:
```python
# PCA figure
showMePCAFig(DataSlice, targetCol)

# Violin plots
showMeViolinPlts(targetCol, signifTVals, DataSlice, boolean, number)

# ROI Accuracy Plot
showMeROIAccPlot(maxROI, Orig_TS_path_names, tsData, targetCol, boolean, indices2Keep)

# Feature Accuracy Plot
showMeFeatAccPlot(featMat3D, targetCol, boolean)

### Extra plots

# The distribution of the fdAvgs
fdAvgs.hist(column='FD Avgs')

# The number of subjects remaining when filtered by decreasing fd
fdArray.plot.line(x='fd', y=['SCZ', 'Control'])

# The variation of balanced accuracy with respect to threshold fd values

```
Too much info
