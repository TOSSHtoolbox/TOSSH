function [results] = calc_Euser(Q_mat, t_mat)
%calc_Euser calculates signatures from Euser et al. (2013).
%   Euser et al. (2013) use 8 signatures that represent different aspects
%   of hydrological behaviour. The signatures are used to test the 
%   consistency of model performance, within the FARM model evaluation 
%   framework.
%   Some of the 8 signatures used in Euser et al. (2013) are the same but
%   applied to different parts of the time series, e.g. the low flow
%   period. There are two ways of doing this. We can either set all values
%   outside the chosen period to NaN, or we can remove them. The latter is
%   the default option here. Since this leads to gaps in the time series 
%   which would raise a warning, we temporarily disable warnings.
%
%   INPUT
%   Q_mat: streamflow [mm/timestep] matrix (cell array)
%   t_mat: time [Matlab datenum] matrix (cell array)
%
%   OUTPUT
%   results: struc array with all results (each signature for each time
%       series and associated error strings)
%
%   EXAMPLE
%   % load example data
%   data = load('example/example_data/33029_daily.mat');
%   % create consistent cell arrays
%   Q_mat = {data.Q};
%   t_mat = {data.t};
%   results = calc_Euser(Q_mat,t_mat);
%
%   References
%   Euser, T., Winsemius, H.C., Hrachowitz, M., Fenicia, F., Uhlenbrook, S.
%   and Savenije, H.H.G., 2013. A framework to assess the realism of model
%   structures using hydrological signatures. Hydrology and Earth System
%   Sciences, 17 (5), 2013.
%
%   Copyright (C) 2020
%   This software is distributed under the GNU Public License Version 3.
%   See <https://www.gnu.org/licenses/gpl-3.0.en.html> for details.

% check input parameters
if nargin < 2
    error('Not enough input arguments.')
end

ip = inputParser;
ip.CaseSensitive = true;

% required input arguments
% Please input time series as a cell array of the following format:
% {x_1; x_2; ...; x_n}, where each entry (1, 2, ..., n) corresponds to one
% time series, e.g. from one catchment. For one catchment only, please
% input {x}. Example: {Q_1; Q_2; ...; Q_n} for streamflow.
addRequired(ip, 'Q_mat', @(Q_mat) iscell(Q_mat))
addRequired(ip, 't_mat', @(t_mat) iscell(t_mat))

parse(ip, Q_mat, t_mat)

% calculate signatures

% initialise arrays
AC1 = NaN(size(Q_mat,1),1);
AC1_error_str = strings(size(Q_mat,1),1);
AC1_low = NaN(size(Q_mat,1),1); % low flow season
AC1_low_error_str = strings(size(Q_mat,1),1);
RLD = NaN(size(Q_mat,1),1);
RLD_error_str = strings(size(Q_mat,1),1);
PeakDistribution = NaN(size(Q_mat,1),1);
PeakDistribution_error_str = strings(size(Q_mat,1),1);
PeakDistribution_low = NaN(size(Q_mat,1),1); % low flow season
PeakDistribution_low_error_str = strings(size(Q_mat,1),1);
FDC = cell(size(Q_mat,1),1); % they use the total FDC to evaluate models
FDC_error_str = strings(size(Q_mat,1),1);
FDC_low = cell(size(Q_mat,1),1); % low flow season
FDC_low_error_str = strings(size(Q_mat,1),1);
FDC_high = cell(size(Q_mat,1),1); % high flow season
FDC_high_error_str = strings(size(Q_mat,1),1);

low_flow_season = [5:9]; % May to September
high_flow_season = [11:12, 1:4]; % November to April

% loop over all catchments
for i = 1:size(Q_mat,1)
    
    [Q_low, t_low] = ...
        util_ExtractSubPeriod(Q_mat{i}, t_mat{i}, low_flow_season);
    [Q_high, t_high] = ...
        util_ExtractSubPeriod(Q_mat{i}, t_mat{i}, high_flow_season);
    
    [AC1(i),~,AC1_error_str(i)] = sig_Autocorrelation(Q_mat{i},t_mat{i},'lag',1);
    [AC1_low(i),~,AC1_low_error_str(i)] = sig_Autocorrelation(Q_low,t_low,'lag',1);
    [RLD(i),~,RLD_error_str(i)] = sig_RisingLimbDensity(Q_mat{i},t_mat{i});
    [PeakDistribution(i),~,PeakDistribution_error_str(i)] = sig_PeakDistribution(Q_mat{i},t_mat{i});
    [PeakDistribution_low(i),~,PeakDistribution_low_error_str(i)] = sig_PeakDistribution(Q_low,t_low);
    [FDC{i}(:,1),FDC{i}(:,2),~,FDC_error_str(i)] = sig_FDC(Q_mat{i},t_mat{i});
    [FDC_low{i}(:,1),FDC_low{i}(:,2),~,FDC_low_error_str(i)] = sig_FDC(Q_low,t_low);
    [FDC_high{i}(:,1),FDC_high{i}(:,2),~,FDC_high_error_str(i)] = sig_FDC(Q_high,t_high);
    
end

% add results to struct array
results.AC1 = AC1;
results.AC1_error_str = AC1_error_str;
results.AC1_low = AC1_low;
results.AC1_low_error_str = AC1_low_error_str;
results.RLD = RLD;
results.RLD_error_str = RLD_error_str;
results.PeakDistribution = PeakDistribution;
results.PeakDistribution_error_str = PeakDistribution_error_str;
results.PeakDistribution_low = PeakDistribution_low;
results.PeakDistribution_low_error_str = PeakDistribution_low_error_str;
results.FDC = FDC;
results.FDC_error_str = FDC_error_str;
results.FDC_low = FDC_low;
results.FDC_low_error_str = FDC_low_error_str;
results.FDC_high = FDC_high;
results.FDC_high_error_str = FDC_high_error_str;

end
