function [results] = calc_BasicSet(Q_mat, t_mat)
%calc_BasicSet calculates basic set of signatures.
%   The basic set of signatures are designed to cover the five components 
%   of a natural streamflow regime as defined by Poff et al. (1997) and 
%   Richter et al. (1996): magnitude, frequency, duration, timing and rate
%   of change. As Poff et al. state, these components "can be used to 
%   characterize the entire range of flows and specific hydrologic 
%   phenomena, such as floods or low flows, that are critical to the 
%   integrity of river ecosystems".
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
%   results = calc_BasicSet(Q_mat,t_mat);
%
%   References
%   Poff, N.L., Allan, J.D., Bain, M.B., Karr, J.R., Prestegaard, K.L.,
%   Richter, B.D., Sparks, R.E. and Stromberg, J.C., 1997. The natural flow
%   regime. BioScience, 47(11), pp.769-784.
%   Richter, B.D., Baumgartner, J.V., Powell, J. and Braun, D.P., 1996. A
%   method for assessing hydrologic alteration within ecosystems.
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
Q_mean = NaN(size(Q_mat,1),1);
Q_mean_error_str = strings(size(Q_mat,1),1);
Q5 = NaN(size(Q_mat,1),1);
Q5_error_str = strings(size(Q_mat,1),1);
Q95 = NaN(size(Q_mat,1),1);
Q95_error_str = strings(size(Q_mat,1),1);
Q_mean_monthly = NaN(size(Q_mat,1),12);
Q_mean_monthly_error_str = strings(size(Q_mat,1),1);
Q_7_day_min = NaN(size(Q_mat,1),1);
Q_7_day_min_error_str = strings(size(Q_mat,1),1);
BFI = NaN(size(Q_mat,1),1);
BFI_error_str = strings(size(Q_mat,1),1);
CoV = NaN(size(Q_mat,1),1);
CoV_error_str = strings(size(Q_mat,1),1);
x_Q_frequency = NaN(size(Q_mat,1),1);
x_Q_frequency_error_str = strings(size(Q_mat,1),1);
x_Q_duration = NaN(size(Q_mat,1),1);
x_Q_duration_error_str = strings(size(Q_mat,1),1);
HFD_mean = NaN(size(Q_mat,1),1);
HFD_mean_error_str = strings(size(Q_mat,1),1);
HFI_mean = NaN(size(Q_mat,1),1);
HFI_mean_error_str = strings(size(Q_mat,1),1);
AC1 = NaN(size(Q_mat,1),1);
AC1_error_str = strings(size(Q_mat,1),1);
FDC_slope = NaN(size(Q_mat,1),1);
FDC_slope_error_str = strings(size(Q_mat,1),1);
BaseflowRecessionK = NaN(size(Q_mat,1),1);
BaseflowRecessionK_error_str = strings(size(Q_mat,1),1);

% loop over all catchments
for i = 1:size(Q_mat,1)
    
    [Q_mean(i),~,Q_mean_error_str(i)] = sig_Q_mean(Q_mat{i},t_mat{i});
    [Q5(i),~,Q5_error_str(i)] = sig_x_percentile(Q_mat{i},t_mat{i},[5]);
    [Q95(i),~,Q95_error_str(i)] = sig_x_percentile(Q_mat{i},t_mat{i},[95]);
    [Q_mean_monthly(i,:),~,Q_mean_monthly_error_str(i)] = sig_Q_mean_monthly(Q_mat{i},t_mat{i},[1:12]);
    [Q_7_day_min(i),~,Q_7_day_min_error_str(i)] = sig_Q_n_day_min(Q_mat{i},t_mat{i},7);
    [BFI(i),~,BFI_error_str(i)] = sig_BFI(Q_mat{i},t_mat{i});
    [CoV(i),~,CoV_error_str(i)] = sig_Q_CoV(Q_mat{i},t_mat{i});
    [x_Q_frequency(i),~,x_Q_frequency_error_str(i)] = sig_x_Q_frequency(Q_mat{i},t_mat{i},'low');
    [x_Q_duration(i),~,x_Q_duration_error_str(i)] = sig_x_Q_duration(Q_mat{i},t_mat{i},'low');
    [HFD_mean(i),~,HFD_mean_error_str(i)] = sig_HFD_mean(Q_mat{i},t_mat{i});
    [HFI_mean(i),~,HFI_mean_error_str(i)] = sig_HFI_mean(Q_mat{i},t_mat{i});
    [AC1(i),~,AC1_error_str(i)] = sig_Autocorrelation(Q_mat{i},t_mat{i},'lag',1);
    [FDC_slope(i),~,FDC_slope_error_str(i)] = sig_FDC_slope(Q_mat{i},t_mat{i});
    [BaseflowRecessionK(i),~,BaseflowRecessionK_error_str(i)] = sig_BaseflowRecessionK(Q_mat{i},t_mat{i},'eps',0.001*median(Q_mat{i},'omitnan'));
    
end

% add results to struct array
results.Q_mean = Q_mean;
results.Q_mean_error_str = Q_mean_error_str;
results.Q5 = Q5;
results.Q5_error_str = Q5_error_str;
results.Q95 = Q95;
results.Q95_error_str = Q95_error_str;
results.Q_mean_monthly = Q_mean_monthly;
results.Q_mean_monthly_error_str = Q_mean_monthly_error_str;
results.Q_7_day_min = Q_7_day_min;
results.Q_7_day_min_error_str = Q_7_day_min_error_str;
results.BFI = BFI;
results.BFI_error_str = BFI_error_str;
results.CoV = CoV;
results.CoV_error_str = CoV_error_str;
results.x_Q_frequency = x_Q_frequency;
results.x_Q_frequency_error_str = x_Q_frequency_error_str;
results.x_Q_duration = x_Q_duration;
results.x_Q_duration_error_str = x_Q_duration_error_str;
results.HFD_mean = HFD_mean;
results.HFD_mean_error_str = HFD_mean_error_str;
results.HFI_mean = HFI_mean;
results.HFI_mean_error_str = HFI_mean_error_str;
results.AC1 = AC1;
results.AC1_error_str = AC1_error_str;
results.FDC_slope = FDC_slope;
results.FDC_slope_error_str = FDC_slope_error_str;
results.BaseflowRecessionK = BaseflowRecessionK;
results.BaseflowRecessionK_error_str = BaseflowRecessionK_error_str;

end
