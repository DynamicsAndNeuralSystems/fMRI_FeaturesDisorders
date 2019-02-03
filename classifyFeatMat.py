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

# Run below code (commented out) in a for loop to generate the fdArray.txt file

df = pd.DataFrame({'fd': [], 'SCZ': [], 'Control': [], 'Total': [],
'SCZ:Control': [], 'AvgROIAcc': [], 'AvgROIError': [], 'AvgFeatAcc': [], 'AvgFeatError': []})

# df = pd.DataFrame({'fd': [], 'feat3BalancedAcc': [], 'stdDev': []})

# Make fd array
filePath = '/Users/AV/Dropbox/UCLA/movementData/fdAvgs.txt'
fdAvgs = pd.read_csv(filePath,header=None);
maxFd =  "%.2f" % fdAvgs.max()
minFd = "%.2f" % fdAvgs.min()
# print(maxFd)
fdArray = np.linspace(0.72,0.12,21)
print(fdArray)

for threshold_fd in fdArray:

    # threshold_fd = 0.20 # Need to remove this if running the code above

    # Store the path of the folder containing the data
    # Alternatively, could use user input -
    # path = input('Enter the file path of the data: ')
    path = '/Users/AV/Dropbox/UCLA/cfgData/'

    # Need to alphabetise and store the filenames in the path
    Orig_TS_path_names = sorted(glob.glob(path + '*.mat'))
    TS_path_names = sorted(glob.glob(path + '*.mat'))

    # In[2]:

    TS_path_names, indices2Keep = af.removePathNames(threshold_fd, TS_path_names)
    indices2Keep = indices2Keep.tolist()

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
    ROISlice, ROI, maxROI = af.getROISlice(Orig_TS_path_names, tsData, 1, indices2Keep)

    # Assign the data to variables
    y = np.ravel(targetCol)

    # In[6]:

    # Feature selection

    # Need to store the data as a 3D feature matrix which will allow access to relevant slices
    # [rows,cols,layers] = i, j, k = 185, 22, ROIs - varies element to element
    featMat3D = np.zeros((len(TS_path_names),22,maxROI))

    # Loop through each slice and store it in the matrix
    for i in range(1, maxROI+1):
        featMat3D[:,:,i-1] = af.getROISlice(Orig_TS_path_names, tsData, i, indices2Keep)[0]

    FeatSlice = pd.DataFrame(af.getFeatSlice(featMat3D,3))

    DataSlice = FeatSlice # or ROISlice
    DataSlice_zscored = DataSlice.apply(zscore)
    X = DataSlice_zscored

    # In[7]:

    # Perform 10-Fold Cross Validation

    # Store the function's output as a variable
    # scores = af.get10FoldCVScore(X,y)

#     feat3BalancedAcc = af.get10FoldCVScore(X,y).mean()
#     stdDev = af.get10FoldCVScore(X,y).std()
#
#     df = df.append({'fd': threshold_fd, 'feat3BalancedAcc': feat3BalancedAcc, 'stdDev': stdDev}, ignore_index=True)
#
#     print(df)
#
# df.to_csv('feat3BalancedAcc.txt', index=False)



    # Print scores
    # print('10-fold CV scores as a percentage: ' + str(scores))
    # print('')

    # Mean 10-fold CV score with an error of 1 std dev
    # print("Accuracy as a percentage: %0.1f (+/- %0.1f)" % (scores.mean(), scores.std()))

    # In[8]:

    # Store the function's main output
    signifTVals = af.getTPVals(targetCol, DataSlice)[2]

    # In[9]:

    # PCA

    # af.showMePCAFig(DataSlice, targetCol)

    # In[10]:

    # Violin plots

    # af.showMeViolinPlts(targetCol, signifTVals, DataSlice, 0, 1)

    # In[11]:

    # Plotting the number of regions with classification accuracies greater than 60%

    AvgROIAcc, AvgROIError = af.showMeROIAccPlot(maxROI, Orig_TS_path_names, tsData, targetCol, 0, indices2Keep)
    # af.showMeROIAccPlot(maxROI, Orig_TS_path_names, tsData, targetCol, 1, indices2Keep)

    print('Avg Reg Acc (%) = ' + str(AvgROIAcc))

    # In[12]:

    # Plot feature accuracies

    AvgFeatAcc, AvgFeatError = af.showMeFeatAccPlot(featMat3D, targetCol, 0)
    # af.showMeFeatAccPlot(featMat3D, targetCol, 1)

    print('Avg Feat Acc (%) = ' + str(AvgFeatAcc))
    print('')

    df = df.append({'fd': threshold_fd, 'SCZ': SCZ, 'Control': Control, 'Total': Total, 'SCZ:Control': SCZ2Ctrl,
                    'AvgROIAcc': AvgROIAcc, 'AvgROIError': AvgROIError, 'AvgFeatAcc': AvgFeatAcc, 'AvgFeatError': AvgFeatError}, ignore_index=True)

    print(df)

df.to_csv('fdArray.txt', index=False)

'''

# In[13]:

#-------------------------------------------------------------------------------
# Plot the fd vs classification accuracy

# Import and store the fdArray
filePath = '/Users/AV/Desktop/FeatureMatrixData/fdArray.txt'
fdArray = pd.read_csv(filePath);

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

plt.xlim(max(fd)+0.02, min(fd)+0.02)
ax1.set_xlabel('fd (mm)')
ax1.set_ylabel('Balanced Acc (%)')

ax2 = ax1.twinx()
ax2.plot(fd, SCZ2Control,'g-')
ax2.set_ylabel('SCZ:Control', color='g')
ax2.tick_params('y', colors='g')
plt.grid(b=None)

ax1.legend(['Avg Feat Acc','Avg ROI Acc'],loc=1)
ax2.legend(['SCZ:Control'],loc=2)
fig.tight_layout()
plt.show()
#-------------------------------------------------------------------------------
# fd variation across subjects

filePath = '/Users/AV/Dropbox/UCLA/movementData/fdAvgs.txt'
fdAvgs = pd.read_csv(filePath,header=None, names=['FD Avgs']);

fdAvgs.hist(column='FD Avgs')

plt.xlim(max(fd), min(fd))
plt.title('FD Distribution')
plt.xlabel('FD Averages (mm)')
plt.ylabel('No. of Subjects')
plt.show()
#-------------------------------------------------------------------------------
# No. of subjects left as fd decreases

SCZ = fdArray.iloc[:,1]
Control = fdArray.iloc[:,2]

fdArray.plot.line(x='fd', y=['SCZ', 'Control'])

plt.xlim(max(fd), min(fd))
plt.title('FD vs Subjects Remaining')
plt.xlabel('FD Averages (mm)')
plt.ylabel('No. of Subjects')
plt.show()
#-------------------------------------------------------------------------------
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
