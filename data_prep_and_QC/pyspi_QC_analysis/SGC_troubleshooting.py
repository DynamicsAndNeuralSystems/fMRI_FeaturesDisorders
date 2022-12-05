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
from pyspi.statistics.spectral import spectral_granger
from pyspi.calculator import Calculator

# Define data path
data_path = "/Users/abry4213/data/UCLA_CNP/raw_data/pydata/"

# Load data that is NaN for the full the frequency spectrum
data_np_fullnan = pd.read_csv(data_path + "sub-10274_lh_lingual_rh_lingual.csv",
                      header=None).to_numpy()
data_np_fullnan_rev = np.flip(data_np_fullnan, 0)

# z-score data
data_np_fullnan_z = scipy.stats.zscore(data_np_fullnan, axis=1)
data_np_fullnan_rev_z = scipy.stats.zscore(data_np_fullnan_rev, axis=1)


###############################################################################
# Case study: sub-10274, NaN for full frequency range
###############################################################################

########### lh_lingual --> rh_lingual ############

# Calculate multitaper and connectivity objects
m_fullnan = Multitaper(
        data_np_fullnan_z.T,
        sampling_frequency=1
)
c_fullnan = Connectivity.from_multitaper(m_fullnan)
    
# Compute pairwise SGC
pairwise_spectral_granger_fullnan=c_fullnan.pairwise_spectral_granger_prediction()

# Obtain positive frequencies
freq_fullnan = m_fullnan.frequencies[m_fullnan.frequencies >= 0]
GA_values_fullnan = c_fullnan.pairwise_spectral_granger_prediction()[0,:,0,1]

# Create dataframe
SGC_df_fullnan = pd.DataFrame({"Freq": freq_fullnan,
                     "SGC": GA_values_fullnan})

# full frequency range mean using spectral-connectivity
np.mean(SGC_df_fullnan.query("Freq > 0")["SGC"].to_numpy())
np.nanmean(SGC_df_fullnan.query("Freq > 0")["SGC"].to_numpy())

# full frequency range mean using pyspi
SGC_pyspi_fullnan = spectral_granger().multivariate(data_np_fullnan_z)
SGC_pyspi_fullnan[0,1]

# lower frequency range mean using spectral-connectivity
np.mean(SGC_df_fullnan.query("Freq <= 0.25 & Freq > 0")["SGC"].to_numpy())
np.nanmean(SGC_df_fullnan.query("Freq <= 0.25 & Freq > 0")["SGC"].to_numpy())

# lower frequency range mean using pyspi
SGC_pyspi_fullnan_lower = spectral_granger(fmin=0, fmax=0.25).multivariate(data_np_fullnan_z)
SGC_pyspi_fullnan_lower[0,1]

# upper frequency range mean using spectral-connectivity
np.mean(SGC_df_fullnan.query("Freq >= 0.25")["SGC"].to_numpy())
np.nanmean(SGC_df_fullnan.query("Freq >= 0.25")["SGC"].to_numpy())

# Upper frequency range mean using pyspi
SGC_pyspi_fullnan_upper = spectral_granger(fmin=0.25, fmax=0.5).multivariate(data_np_fullnan_z)
SGC_pyspi_fullnan_upper[0,1]

# Find NaN values at individual frequencies from spectral-connectivity
SGC_df_fullnan.query("is.nan(SGC)")

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
