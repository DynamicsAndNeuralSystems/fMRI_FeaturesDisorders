function out = SP_fALFF(y,samplingPeriod)
% SP_fALFF  Computes the fALFF statistic popular in fMRI analysis
%
% Implementation based on that in the REST toolkit:
% cf. http://restfmri.net/forum/index.php
%
% Song et al. (2011) REST: A Toolkit for Resting-State Functional Magnetic
% Resonance Imaging Data Processing. PLoS ONE 6(9): e25031.
% doi:10.1371/journal.pone.0025031
%
%---INPUTS:
%
% y, the input time series
%
% samplingPeriod, the sampling period.

% ------------------------------------------------------------------------------
% Copyright (C) 2016, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2013). DOI: 10.1098/rsif.2013.0048
%
% This function is free software: you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation, either version 3 of the License, or (at your option) any later
% version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License along with
% this program. If not, see <http://www.gnu.org/licenses/>.
% ------------------------------------------------------------------------------

if nargin < 2
    samplingPeriod = 2; % (s)
end

%-------------------------------------------------------------------------------
% Parameters:
%-------------------------------------------------------------------------------
% Band pass parameters (default window for ALFF):
lowCutoff = 0.01; % (Hz)
highCutoff = 0.08; % (Hz)

%-------------------------------------------------------------------------------
% Get the frequency index
sampleFreq 	 = 1/samplingPeriod;
sampleLength = length(y);
paddedLength = 2^nextpow2(sampleLength);

freq = 0:(sampleFreq/paddedLength):(sampleFreq/2-sampleFreq/paddedLength);
if lowCutoff >= sampleFreq/2 % All high included
    idx_lowCutoff = paddedLength/2 + 1;
else % high cut off, such as freq > 0.01 Hz
    idx_lowCutoff = ceil(lowCutoff * paddedLength * samplingPeriod + 1);
end
if (highCutoff >= sampleFreq/2) || (highCutoff==0) % All low pass
    idx_highCutoff = paddedLength/2 + 1;
else % Low pass, such as freq < 0.08 Hz
    idx_highCutoff = floor(highCutoff * paddedLength * samplingPeriod + 1);
end

% Detrend before fft as did in the previous f_alff.m
y = detrend(y);

% Zero padding
y = [y; zeros(paddedLength - sampleLength,1)];

% Perform FFT:
Y = (2*abs(fft(y)).^2)/sampleLength;

% % Compute fALFF measure:
out.fALFF = sum(Y(idx_lowCutoff:idx_highCutoff)) / sum(Y(2:(paddedLength/2 + 1)));
%
% % Raw ALFF?:
out.ALFF = mean(Y(idx_lowCutoff:idx_highCutoff));

% Plot:
freqVector = 0:(sampleFreq/paddedLength):(sampleFreq/2-sampleFreq/paddedLength);
oneSidedPower = Y(1:end/2);


%%
% Power ALFF
out.ALFFpower = (sum(2*oneSidedPower).*(sampleFreq/paddedLength));
end
