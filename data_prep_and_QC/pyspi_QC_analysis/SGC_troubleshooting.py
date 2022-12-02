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

# time-series: paired tuples
data_path = "/Users/abry4213/data/UCLA_Schizophrenia/raw_data/pydata/"
data_np = pd.read_csv(data_path + "sub-10527_lh_rostralanteriorcingulate_lh_caudalmiddlefrontal.csv",
                      header=None).to_numpy()
data_np_rev = np.flip(data_np, 0)

# z-score data
data_np_z = scipy.stats.zscore(data_np, axis=1)
data_np_rev_z = scipy.stats.zscore(data_np_rev, axis=1)

###############################################################################
# Non-parametric SGC from spectral_connectivity
###############################################################################

########### lh_rostralanteriorcingulate --> lh_caudalmiddlefrontal ############


# Calculate multitaper and connectivity objects
m = Multitaper(
        data_np_z.T,
        sampling_frequency=1
)
c = Connectivity.from_multitaper(m)
    
# Compute pairwise SGC
pairwise_spectral_granger=c.pairwise_spectral_granger_prediction()

# Obtain positive frequencies
freq = m.frequencies[m.frequencies >= 0]
GA_values = c.pairwise_spectral_granger_prediction()[0,:,0,1]

# Create dataframe
SGC_df = pd.DataFrame({"Freq": freq,
                     "SGC": GA_values})

# full frequency range mean
np.mean(SGC_df.SGC)

# lower frequency range mean
SGC_df.query("Freq <= 0.25")["SGC"].mean()

# Upper frequency range mean
SGC_df.query("Freq >= 0.25")["SGC"].mean()


########### lh_caudalmiddlefrontal --> lh_rostralanteriorcingulate ############

# Calculate multitaper and connectivity objects
m_rev = Multitaper(
        data_np_rev_z.T,
        sampling_frequency=1
)
c_rev = Connectivity.from_multitaper(m_rev)
    
# Compute pairwise SGC
GA_values_rev = c_rev.pairwise_spectral_granger_prediction()[0,:,0,1]

# Obtain positive frequencies
freq_rev = m_rev.frequencies[m_rev.frequencies >= 0]

# Create dataframe
SGC_df_rev = pd.DataFrame({"Freq": freq_rev,
                     "SGC": GA_values_rev})

# full frequency range mean
np.mean(SGC_df_rev.SGC)

# lower frequency range mean
SGC_df_rev.query("Freq <= 0.25")["SGC"].mean()

# Upper frequency range mean
SGC_df_rev.query("Freq >= 0.25")["SGC"].mean()

###############################################################################
# Non-parametric SGC directly from pyspi
###############################################################################

from pyspi.statistics.spectral import spectral_granger
from pyspi.calculator import Calculator

# full frequency range 
SGC_pyspi = spectral_granger().multivariate(data_np_z)

# Lower frequency range
SGC_pyspi_lower = spectral_granger(fmin=0, fmax=0.25).multivariate(data_np_z)

# Upper frequency range
SGC_pyspi_upper = spectral_granger(fmin=0.25, fmax=0.5).multivariate(data_np_z)

###############################################################################
# Parametric SGC from nitime
###############################################################################

# Instantiate time-series
time_series = ts.TimeSeries(data_np_z,sampling_interval=1)

# ## parametric
# GA = nta.GrangerAnalyzer(time_series, order=1)
# GA_values = GA.causality_xy[0, 1, :]
# GA_values_rev = GA.causality_yx[0, 1, :]
# freq = GA.frequencies
