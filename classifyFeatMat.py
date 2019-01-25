#!/usr/local/bin/python3

# In[1]:

# Setting up

# Import required modules
import glob
import scipy.io as sio
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('ggplot')
import seaborn as sns
import sklearn
import analysisFunctions as af

'''
Run below code (commented out) in a for loop to generate the fdArray.txt file

# df = pd.DataFrame({'fd': [], 'Total': [], 'SCZ:Control': [], 'AvgRegAcc': [], 'AvgFeatAcc': []})
#
# # Make fd array
# fdArray = np.linspace(0.2,0.12,9)
# print(fdArray)
#
# for threshold_fd in fdArray:
'''

threshold_fd = 0.16 # Need to remove this if running the code above

# Store the path of the folder containing the data
# Alternatively, could use user input -
# path = input('Enter the file path of the data: ')
path = '/Users/AV/Dropbox/UCLA/cfgData/'

# Need to alphabetise and store the filenames in the path
TS_path_names = sorted(glob.glob(path + '*.mat'))

# In[2]:

TS_path_names = af.removePathNames(threshold_fd, TS_path_names)

# In[3]:

# Reading in the feature matrix data from the .txt file

# Import and store the feature matrix in the variable tsData
filePath = '/Users/AV/Desktop/FeatureMatrixData/element1.txt'
tsData = pd.read_csv(filePath,header=None);

# In[4]:

# Creating the target column

targetCol = af.getTargetCol(TS_path_names)

Control = (targetCol == 0).sum()
print('Control = ' + str(Control))
SCZ = (targetCol == 1).sum()
print('SCZ = ' + str(SCZ))
print('')
Total = int(SCZ + Control)
SCZ2Ctrl = '{0:.2f}'.format(SCZ/Control)

# In[5]:

# ROI selection

# Need to z-score the selection of tsData
from scipy.stats import zscore
tsDataSlice, ROI, maxROI = af.getROISlice(TS_path_names, tsData, 1)
# tsDataSlice_zscored = tsDataSlice.apply(zscore)

# # Assign the data to variables
# X = tsDataSlice_zscored
# y = np.ravel(targetCol)

# In[6]:

# Need to store the data as a 3D feature matrix which will allow access to relevant slices
# [rows,cols,layers] = i, j, k = 185, 22, ROIs - varies element to element

featMat3D = np.zeros((len(TS_path_names),22,maxROI))

# Loop through each slice and store it in the matrix
for i in range(1, maxROI+1):
    featMat3D[:,:,i-1] = af.getROISlice(TS_path_names, tsData, i)[0]

# In[7]:

# Perform 10-Fold Cross Validation

# Store the function's output as a variable
# scores = af.get10FoldCVScore(X,y)

# # Print scores
# print('10-fold CV scores as a percentage: ' + str(scores))
# print('')
#
# # Mean 10-fold CV score with an error of 1 std dev
# print("Accuracy as a percentage: %0.1f (+/- %0.1f)" % (scores.mean(), scores.std()))

# In[8]:

# Store the function's main output
signifTVals = af.getTPVals(targetCol, tsDataSlice)[2]

# In[9]:

# PCA

# af.showMePCAFig(tsDataSlice, targetCol)

# In[10]:

# Violin plots

# af.showMeViolinPlts(targetCol, signifTVals, X, ROI)

# In[11]:

# Plotting the number of regions with classification accuracies greater than 60%

''' # AvgRegAcc = '''
af.showMeRegAccPlot(maxROI, TS_path_names, tsData, targetCol, 1)
# print('Avg Reg Acc (%) = ' + str(AvgRegAcc))

# In[12]:

# Plot feature accuracies

''' # AvgFeatAcc = '''
af.showMeFeatAccPlot(featMat3D, targetCol, 1)

'''
# print('Avg Feat Acc (%) = ' + str(AvgFeatAcc))
# print('')

# df = df.append({'fd': threshold_fd, 'Total': Total, 'SCZ:Control': SCZ2Ctrl,
#                 'AvgRegAcc': AvgRegAcc, 'AvgFeatAcc': AvgFeatAcc}, ignore_index=True)

# print(df)

# df.to_csv('fdArray.txt')
'''

'''
For reference ONLY:

List of catch22 features
1.	CO_Embed2_Dist_tau_d_expfit_meandiff
2.	CO_FirstMin_ac
3.	CO_HistogramAMI_even_2_5
4.	CO_f1ecac
5.	CO_trev_1_num
6.	DN_HistogramMode_10
7.	DN_HistogramMode_5
8.	DN_OutlierInclude_n_001_mdrmd
9.	DN_OutlierInclude_p_001_mdrmd
10.	FC_LocalSimple_mean1_tauresrat
11.	FC_LocalSimple_mean3_stderr
12.	IN_AutoMutualInfoStats_40_gaussian_fmmi
13.	MD_hrv_classic_pnn40
14.	PeriodicityWang_th0_01
15.	SB_BinaryStats_diff_longstretch0
16.	SB_BinaryStats_mean_longstretch1
17.	SB_MotifThree_quantile_hh
18.	SB_TransitionMatrix_3ac_sumdiagcov
19.	SC_FluctAnal_2_dfa_50_1_2_logi_prop_r1
20.	SC_FluctAnal_2_rsrangefit_50_1_logi_prop_r1
21.	SP_Summaries_welch_rect_area_5_1
22.	SP_Summaries_welch_rect_centroid
'''
