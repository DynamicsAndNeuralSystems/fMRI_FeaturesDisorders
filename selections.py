import glob
import after_catch_analysis_pack as acap
import pandas as pd
import numpy as np
from numpy import genfromtxt

'''Use this file to make all your selections e.g. which dataset will be used,
paths to files etc.'''

class Selections:
    def __init__(self, dataset=""):
        #------------------------MAKE SELECTIONS HERE---------------------------
        if len(dataset) == 0:
            self.dataset = "UCLA" # Select which dataset is being used
        else:
            self.dataset=dataset
        c22DataFileName = 'c22_roiTS1_'+self.dataset+'.txt'

        # Select path to folder containing the subject data for the given dataset
        folderPath = '/Users/preethompal/Dropbox (Sydney Uni Student)/'+self.dataset+'/cfgData/'

        # Select path to framewise displacement data.
        fdAvgsPath = '/Users/preethompal/Dropbox (Sydney Uni Student)/'+self.dataset+'/movementData/fdAvgs_'+self.dataset+'.txt'

        # Select path to participants.csv file.
        csvPath = '/Users/preethompal/Dropbox (Sydney Uni Student)/'+self.dataset+'/participants.csv'

        featureListFileName = 'PythonFeatureList.txt' # This text file contains the names of the 22 features.

        # Choose whether or not to display plots for each feature and region while running.
        # Will slow down but may be interesting.
        self.dispFigs = False
        #-----------------------------------------------------------------------

        #-----------------------DON'T CHANGE BELOW HERE-------------------------
        self.featNames = [lines.rstrip('\n') for lines in open(featureListFileName)]
        self.participants = pd.read_csv(csvPath,header=0);
        self.filePathsAll = sorted(glob.glob(folderPath + '*.mat'))
        self.subjCount = len(self.filePathsAll) #total number of subjects
        c22DataNdArray = genfromtxt(c22DataFileName, delimiter=',')
        # Convert c22Data to a data frame with a MultiIndex. Get total number of regions and features
        self.c22Data, self.roiCount, self.featCount = acap.addMultiIndex(c22DataNdArray,self.featNames, self.subjCount)
        self.fdAvgs = pd.read_csv(fdAvgsPath,header=None, names=['Avgs']);
