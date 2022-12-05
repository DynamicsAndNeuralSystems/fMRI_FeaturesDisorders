#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  1 23:36:09 2022

@author: abry4213
"""

import nitime.analysis as nta
import pandas as pd
import numpy as np
import nitime.timeseries as ts
import scipy.stats
from spectral_connectivity import Multitaper, Connectivity
import warnings

# Define data path
data_path = "/Users/abry4213/data/UCLA_Schizophrenia/raw_data/pydata/"

# Load data that is NaN at the low end of the frequency spectrum
data_np_lownan = pd.read_csv(data_path + "sub-10527_lh_rostralanteriorcingulate_lh_caudalmiddlefrontal.csv",
                      header=None).to_numpy()
data_np_lownan_rev = np.flip(data_np_lownan, 0)

# z-score data
data_np_lownan_z = scipy.stats.zscore(data_np_lownan, axis=1)
data_np_lownan_rev_z = scipy.stats.zscore(data_np_lownan_rev, axis=1)

# Load data that is NaN at the high end of the frequency spectrum
data_np_highnan = pd.read_csv(data_path + "sub-10206_rh_precentral_rh_postcentral.csv",
                      header=None).to_numpy()
data_np_highnan_rev = np.flip(data_np_highnan, 0)

# z-score data
data_np_highnan_z = scipy.stats.zscore(data_np_highnan, axis=1)
data_np_highnan_rev_z = scipy.stats.zscore(data_np_highnan_rev, axis=1)


###############################################################################
# Case study: sub-10527, NaN for low frequnecy range
###############################################################################

########### lh_rostralanteriorcingulate --> lh_caudalmiddlefrontal ############


# Calculate multitaper and connectivity objects
m_lownan = Multitaper(
        data_np_lownan_z.T,
        sampling_frequency=1
)
c_lownan = Connectivity.from_multitaper(m_lownan)
    
# Compute pairwise SGC
pairwise_spectral_granger_lownan=c_lownan.pairwise_spectral_granger_prediction()

# Obtain positive frequencies
freq_lownan = m_lownan.frequencies[m_lownan.frequencies >= 0]
GA_values_lownan = c_lownan.pairwise_spectral_granger_prediction()[0,:,0,1]

# Create dataframe
SGC_df_lownan = pd.DataFrame({"Freq": freq_lownan,
                     "SGC": GA_values_lownan})

# full frequency range mean
np.mean(SGC_df_lownan.SGC)

# lower frequency range mean
SGC_df_lownan.query("Freq <= 0.25")["SGC"].mean()

# Upper frequency range mean
SGC_df_lownan.query("Freq >= 0.25")["SGC"].mean()


########### lh_caudalmiddlefrontal --> lh_rostralanteriorcingulate ############

# Calculate multitaper and connectivity objects
m_lownan_rev = Multitaper(
        data_np_lownan_rev_z.T,
        sampling_frequency=1
)
c_lownan_rev = Connectivity.from_multitaper(m_lownan_rev)
    
# Compute pairwise SGC
GA_values_lownan_rev = c_lownan_rev.pairwise_spectral_granger_prediction()[0,:,0,1]

# Obtain positive frequencies
freq_lownan_rev = m_lownan_rev.frequencies[m_lownan_rev.frequencies >= 0]

# Create dataframe
SGC_df_lownan_rev = pd.DataFrame({"Freq": freq_lownan_rev,
                     "SGC": GA_values_lownan_rev})

# full frequency range mean
np.mean(SGC_df_lownan_rev.SGC)

# lower frequency range mean
SGC_df_lownan_rev.query("Freq <= 0.25")["SGC"].mean()

# Upper frequency range mean
SGC_df_lownan_rev.query("Freq >= 0.25")["SGC"].mean()


###### Non-parametric SGC directly from pyspi
from pyspi.statistics.spectral import spectral_granger
from pyspi.calculator import Calculator

# full frequency range 
SGC_pyspi_lownan = spectral_granger().multivariate(data_np_lownan_z)

# Lower frequency range
SGC_pyspi_lownan_lower = spectral_granger(fmin=0, fmax=0.25).multivariate(data_np_lownan_z)

# Upper frequency range
SGC_pyspi_lownan_upper = spectral_granger(fmin=0.25, fmax=0.5).multivariate(data_np_lownan_z)



###############################################################################
# Case study: sub-10527, NaN for low frequnecy range
###############################################################################

########### lh_rostralanteriorcingulate --> lh_caudalmiddlefrontal ############

# Calculate multitaper and connectivity objects
m_highnan = Multitaper(
        data_np_highnan_z.T,
        sampling_frequency=1
)
c_highnan = Connectivity.from_multitaper(m_highnan)
    
# Compute pairwise SGC
pairwise_spectral_granger_highnan=c_highnan.pairwise_spectral_granger_prediction()

# Obtain positive frequencies
freq_highnan = m_highnan.frequencies[m_highnan.frequencies >= 0]
GA_values_highnan = c_highnan.pairwise_spectral_granger_prediction()[0,:,0,1]

# Create dataframe
SGC_df_highnan = pd.DataFrame({"Freq": freq_highnan,
                     "SGC": GA_values_highnan})

# full frequency range mean
np.mean(SGC_df_highnan.SGC)

# lower frequency range mean
SGC_df_highnan.query("Freq <= 0.25")["SGC"].mean()

# Upper frequency range mean
SGC_df_highnan.query("Freq >= 0.25")["SGC"].mean()


###### Non-parametric SGC directly from pyspi
from pyspi.statistics.spectral import spectral_granger
from pyspi.calculator import Calculator

# full frequency range 
SGC_pyspi_highnan = spectral_granger().multivariate(data_np_highnan_z)

# Lower frequency range
SGC_pyspi_highnan_lower = spectral_granger(fmin=0, fmax=0.25).multivariate(data_np_highnan_z)

# Upper frequency range
SGC_pyspi_highnan_upper = spectral_granger(fmin=0.25, fmax=0.5).multivariate(data_np_highnan_z)



###############################################################################
# Parametric SGC from nitime
###############################################################################

# Instantiate time-series
time_series = ts.TimeSeries(data_np_lownan_z,sampling_interval=1)

# ## parametric
# GA = nta.GrangerAnalyzer(time_series, order=1)
# GA_values = GA.causality_xy[0, 1, :]
# GA_values_rev = GA.causality_yx[0, 1, :]
# freq = GA.frequencies
