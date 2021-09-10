function [results] = calc_McMillan_OverlandFlow(Q_mat, t_mat, P_mat, varargin)
%calc_McMillan_OverlandFlow calculates various overland flow signatures.
%   Calculates 10 overland flow (infiltration excess and saturation excess)
%   signatures from McMillan (2020). These signatures come from previous
%   experimental studies that link catchment or hillslope processes to
%   streamflow response dynamics. Some signatures are implemented direct
%   from the original papers, others are interpreted from a qualitative
%   description in the paper.
%
%   INPUT
%   Q_mat: streamflow [mm/timestep] matrix (cell array)
%   t_mat: time [Matlab datenum] matrix (cell array)
%   P_mat: precipitation [mm/timestep] matrix (cell array)
%   OPTIONAL
%   plot_results: whether to plot results, default = false
%   max_recessiondays: maximum length of recession period after rainfall
%       [days] deemed to be part of the preceeding event (shorter if
%       another event starts), default = 1
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
%   results = calc_McMillan_OverlandFlow(Q_mat,t_mat,P_mat);
%
%   References
%   McMillan, H., 2020. Linking hydrologic signatures to hydrologic
%   processes: A review. Hydrological Processes, 34(6), pp.1393-1409.
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

% optional input arguments
addParameter(ip, 'plot_results', false, @islogical) % whether to plot results
addParameter(ip, 'max_recessiondays', 1, @isnumeric) % maximum number of days to allow recession after rain

parse(ip, Q_mat, t_mat, P_mat, varargin{:})
plot_results = ip.Results.plot_results;
max_recessiondays = ip.Results.max_recessiondays;

% initialise arrays

% Infiltration and saturation excess importance, based on their average
% coefficients in regression equations to predict event flow
% characteristics, adapted from qualitative description in Estrany et al.
% (2010).
IE_effect = NaN(size(Q_mat,1),1);
SE_effect = NaN(size(Q_mat,1),1);

% Significance (using likelihood ratio test) and location of a threshold in
% a plot of quickflow volume vs. maximum intensity, signifying IE process
% (Ali et al., 2013). IE is indicated when IE_thresh_sig < 0.05.
IE_thresh_signif = NaN(size(Q_mat,1),1);
IE_thresh = NaN(size(Q_mat,1),1);

% Significance, location and above-threshold slope of a threshold in a plot
% of quickflow volume vs. total precipitation. SE is indicated when
% SE_thresh_sig < 0.05. Where there is no threshold, indicates flow
% generation from riparian areas (Tani, 1997). Slope above threshold
% indicates rate at which saturated areas expand (Tani, 1997; Becker and
% McDonnell 1998).
SE_thresh_signif = NaN(size(Q_mat,1),1);
SE_thresh = NaN(size(Q_mat,1),1);
SE_slope = NaN(size(Q_mat,1),1);

% Significance and location of a threshold in a plot of quickflow volume
% vs. antecedent precipitation index + total precipitation. SE is indicated
% when storage_thresh_sig < 0.05 (Ali et al., 2013; McGrath et al., 2007).
Storage_thresh_signif = NaN(size(Q_mat,1),1);
Storage_thresh = NaN(size(Q_mat,1),1);

% Minimum quickflow as a percentage of precipitation indicates impermeable
% area contribution (Becker and McDonnell, 1998).
min_Qf_perc = NaN(size(Q_mat,1),1);

% variable to store error strings
OF_error_str = strings(size(Q_mat,1),1);

% loop over all catchments
for i = 1:size(Q_mat,1)
    
    [IE_effect(i),SE_effect(i),IE_thresh_signif(i),IE_thresh(i), ...
        SE_thresh_signif(i),SE_thresh(i),SE_slope(i),Storage_thresh(i), ...
        Storage_thresh_signif(i),min_Qf_perc(i),~,OF_error_str(i)] ...
        = sig_EventGraphThresholds(Q_mat{i},t_mat{i},P_mat{i},...
        'plot_results',plot_results,'max_recessiondays',max_recessiondays);
    
end

% add results to struct array
results.IE_effect = IE_effect;
results.SE_effect = SE_effect;
results.IE_thresh_signif = IE_thresh_signif;
results.SE_thresh_signif = SE_thresh_signif;
results.IE_thresh = IE_thresh;
results.SE_thresh = SE_thresh;
results.SE_slope = SE_slope;
results.Storage_thresh_signif = Storage_thresh_signif;
results.Storage_thresh = Storage_thresh;
results.min_Qf_perc = min_Qf_perc;
results.OF_error_str = OF_error_str;

end
