function [results] = calc_Addor(Q_mat, t_mat, P_mat)
%calc_Addor calculates signatures from Addor et al. (2018).
%   Addor et al. (2018) use 15 signatures that "characterize different
%   parts of the hydrograph, and [...] are sensitive to processes occurring
%   over different time scales". The signatures were selected from those
%   commonly used in the literature, and are used to explore the strength
%   of relationships between signatures and catchment attributes.
%
%   INPUT
%   Q_mat: streamflow [mm/timestep] matrix (cell array)
%   t_mat: time [Matlab datenum] matrix (cell array)
%   P_mat: precipitation [mm/timestep] matrix (cell array)
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
%   P_mat = {data.P}; 
%   results = calc_Addor(Q_mat,t_mat,P_mat);
%
%   References
%   Addor, N., Nearing, G., Prieto, C., Newman, A.J., Le Vine, N. and
%   Clark, M.P., 2018. A ranking of hydrological signatures based on their
%   predictability in space. Water Resources Research, 54(11), pp.8792-8812.
%
%   Copyright (C) 2020
%   This software is distributed under the GNU Public License Version 3.
%   See <https://www.gnu.org/licenses/gpl-3.0.en.html> for details.

% check input parameters
if nargin < 3
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
addRequired(ip, 'P_mat', @(P_mat) iscell(P_mat))

parse(ip, Q_mat, t_mat, P_mat)

% calculate signatures

% initialise arrays
Q_mean = NaN(size(Q_mat,1),1);
Q_mean_error_str = strings(size(Q_mat,1),1);
TotalRR = NaN(size(Q_mat,1),1);
TotalRR_error_str = strings(size(Q_mat,1),1);
QP_elasticity = NaN(size(Q_mat,1),1);
QP_elasticity_error_str = strings(size(Q_mat,1),1);
FDC_slope = NaN(size(Q_mat,1),1);
FDC_slope_error_str = strings(size(Q_mat,1),1);
BFI = NaN(size(Q_mat,1),1);
BFI_error_str = strings(size(Q_mat,1),1);
HFD_mean = NaN(size(Q_mat,1),1);
HFD_mean_error_str = strings(size(Q_mat,1),1);
Q5 = NaN(size(Q_mat,1),1);
Q5_error_str = strings(size(Q_mat,1),1);
Q95 = NaN(size(Q_mat,1),1);
Q95_error_str = strings(size(Q_mat,1),1);
high_Q_freq = NaN(size(Q_mat,1),1);
high_Q_freq_error_str = strings(size(Q_mat,1),1);
high_Q_dur = NaN(size(Q_mat,1),1);
high_Q_dur_error_str = strings(size(Q_mat,1),1);
low_Q_freq = NaN(size(Q_mat,1),1);
low_Q_freq_error_str = strings(size(Q_mat,1),1);
low_Q_dur = NaN(size(Q_mat,1),1);
low_Q_dur_error_str = strings(size(Q_mat,1),1);
zero_Q_freq = NaN(size(Q_mat,1),1);
zero_Q_freq_error_str = strings(size(Q_mat,1),1);

% loop over all catchments
for i = 1:size(Q_mat,1)
    
    [Q_mean(i),~,Q_mean_error_str(i)] = sig_Q_mean(Q_mat{i},t_mat{i});
    [TotalRR(i),~,TotalRR_error_str(i)] = sig_TotalRR(Q_mat{i},t_mat{i},P_mat{i});
    [QP_elasticity(i),~,QP_elasticity_error_str(i)] = ...
        sig_QP_elasticity(Q_mat{i},t_mat{i},P_mat{i},'method','Sanka','start_water_year',10); %,'start_water_year',4 in Southern Hemisphere
    [FDC_slope(i),~,FDC_slope_error_str(i)] = sig_FDC_slope(Q_mat{i},t_mat{i});
    [BFI(i),~,BFI_error_str(i)] = sig_BFI(Q_mat{i},t_mat{i},'method','Lyne_Hollick','parameters',[0.925 3]);
    [HFD_mean(i),~,HFD_mean_error_str(i)] = sig_HFD_mean(Q_mat{i},t_mat{i},'start_month',10); %,'start_month',4 in Southern Hemisphere
    [Q5(i),~,Q5_error_str(i)] = sig_x_percentile(Q_mat{i},t_mat{i},5);
    [Q95(i),~,Q95_error_str(i)] = sig_x_percentile(Q_mat{i},t_mat{i},95);
    [high_Q_freq(i),~,high_Q_freq_error_str(i)] = sig_x_Q_frequency(Q_mat{i},t_mat{i},'high');
    [high_Q_dur(i),~,high_Q_dur_error_str(i)] = sig_x_Q_duration(Q_mat{i},t_mat{i},'high');
    [low_Q_freq(i),~,low_Q_freq_error_str(i)] = sig_x_Q_frequency(Q_mat{i},t_mat{i},'low');
    [low_Q_dur(i),~,low_Q_dur_error_str(i)] = sig_x_Q_duration(Q_mat{i},t_mat{i},'low');
    [zero_Q_freq(i),~, zero_Q_freq_error_str(i)] = sig_x_Q_frequency(Q_mat{i},t_mat{i},'no');
    
end

% add results to struct array
results.Q_mean = Q_mean;
results.Q_mean_error_str = Q_mean_error_str;
results.TotalRR = TotalRR;
results.TotalRR_error_str = TotalRR_error_str;
results.QP_elasticity = QP_elasticity;
results.QP_elasticity_error_str = QP_elasticity_error_str;
results.FDC_slope = FDC_slope;
results.FDC_slope_error_str = FDC_slope_error_str;
results.BFI = BFI;
results.BFI_error_str = BFI_error_str;
results.HFD_mean = HFD_mean;
results.HFD_mean_error_str = HFD_mean_error_str;
results.Q5 = Q5;
results.Q5_error_str = Q5_error_str;
results.Q95 = Q95;
results.Q95_error_str = Q95_error_str;
results.high_Q_freq = high_Q_freq;
results.high_Q_freq_error_str = high_Q_freq_error_str;
results.high_Q_dur = high_Q_dur;
results.high_Q_dur_error_str = high_Q_dur_error_str;
results.low_Q_freq = low_Q_freq;
results.low_Q_freq_error_str = low_Q_freq_error_str;
results.low_Q_dur = low_Q_dur;
results.low_Q_dur_error_str = low_Q_dur_error_str;
results.zero_Q_freq = zero_Q_freq;
results.zero_Q_freq_error_str = zero_Q_freq_error_str;

end
