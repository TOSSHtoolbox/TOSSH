function [Recession_a_Seasonality, error_flag, error_str, fig_handles] ...
    = sig_SeasonalVarRecessions(Q, t, varargin)
%sig_SeasonalVarRecessions calculates seasonal variation in recession parameters.
%   Seasonal variations in recession rate provides information about the
%   impact of evapotranspiration on watershed storage.
%   
%   Notes:
%   Signature only recommended in watersheds with primarily deciduous
%   vegetation in which ET strongly varies seasonally.
%   Original paper uses daily data and this gives more robust results than
%   hourly in testing.
%   Assumes that all individual recession have a slope of 2, and then
%   plots the y-intercept for all individual recession events.
%   Influence of ET in summer is related to the change in intercept with
%   time of year (no exact relationship given).
%
%   INPUT
%   Q: streamflow [mm/timestep]
%   t: time [Matlab datetime]
%   OPTIONAL
%   recession_length: min. length of recession segments [days], default = 5
%   n_start: days to be removed after start of recession
%   eps: allowed increase in flow during recession period, default = 0
%   start_of_recession: define start of recession when baseflow filter
%       rejoins the curve ("baseflow"), or after hydrograph peak ("peak")
%   filter_par: smoothing parameter of Lyne-Hollick filter to determine
%      start of recession (higher = later recession start), default = 0.925
%   plot_results: whether to plot results, default = false
%
%   OUTPUT
%   Recession_a_Seasonality: seasonal change in recession alpha [-]
%       (y-intercept in recession plot assuming a slope of 2)
%   recession_month: approx. month of recession
%   error_flag: 0 (no error), 1 (warning), 2 (error in data check), 3
%       (error in signature calculation)
%   error_str: string contraining error description
%   fig_handles: figure handles to manipulate figures (empty if plotting is
%       not requested)
%
%   EXAMPLE
%   % load example data
%   data = load('example/example_data/33029_daily.mat');
%   Q = data.Q;
%   t = data.t;
%   Recession_a_Seasonality = sig_SeasonalVarRecessions(Q, t);
%   Recession_a_Seasonality = sig_SeasonalVarRecessions(Q, t, 'plot_results', true);
%
%   References
%   Shaw, S. B., & Riha, S. J. (2012). Examining individual recession
%   events instead of a data cloud: Using a modified interpretation
%   of dQ/dt–Q streamflow recession in glaciated watersheds to better
%   inform models of low flow. Journal of Hydrology, 434, 46–54.
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
% time series have to be numeric and either a (n,1) or a (1,n) vector
addRequired(ip, 'Q', @(Q) isnumeric(Q) && (size(Q,1)==1 || size(Q,2)==1))
% date time series has to be numeric or datetime and either a (n,1) or a (1,n) vector
addRequired(ip, 't', @(t) (isnumeric(t) || isdatetime(t)) && (size(t,1)==1 || size(t,2)==1))

% optional input arguments
addParameter(ip, 'recession_length', 5, @isnumeric) % length of decreasing
% flow section (amount of timesteps) to be declared a recession
addParameter(ip, 'n_start', 1, @isnumeric) % days to be removed at beginning of recession
addParameter(ip, 'eps', 0, @isnumeric) % allowed increase in flow during recession period
addParameter(ip, 'start_of_recession', 'peak', @ischar) % defines start of a recession
addParameter(ip, 'filter_par', 0.925, @isnumeric) % smoothing parameter of
% Lyne-Hollick Filter to determine start of recession (higher = later recession start)
addParameter(ip, 'plot_results', false, @islogical) % whether to plot results (2 graphs)
% addParameter(ip, 'fitting_type', 'linear', @ischar) % nonlinear or linear fit

parse(ip, Q, t, varargin{:})
recession_length = ip.Results.recession_length;
n_start = ip.Results.n_start;
eps = ip.Results.eps;
start_of_recession = ip.Results.start_of_recession;
filter_par = ip.Results.filter_par;
plot_results = ip.Results.plot_results;
% fitting_type = ip.Results.fitting_type;

% create empty figure handle
fig_handles = [];

% data checks
[error_flag, error_str, timestep, t] = util_DataCheck(Q, t);
if error_flag == 2
    Recession_a_Seasonality = NaN;
    return
end

% calculate signature

% run recession analysis for individual recessions
error_flag_tmp = 0; % reset error flag since data check will be performed again
error_str_tmp = '';
[para_mat, recession_month, error_flag, error_str, fig_handles] = ...
    sig_RecessionAnalysis(Q, t, ...
    'recession_length', recession_length, ...
    'n_start', n_start, 'eps', eps, 'start_of_recession', start_of_recession, ...
    'filter_par', filter_par, 'fit_individual', true, 'fitting_type', 'slope2');
if error_flag == 3
    Recession_a_Seasonality = NaN;
    return
else
    error_flag = max([error_flag_tmp, error_flag]);
    error_str = [error_str_tmp, error_str];
end

month_median = grpstats(log(para_mat(:,1)),recession_month,{'median'});

Recession_a_Seasonality = max(month_median) - min(month_median);

% optional plotting
if plot_results
    fig = figure('Position',[100 100 350 300]); hold on
    boxplot(log(para_mat(:,1)),recession_month)
    xlabel('Recession Month')
    ylabel('Fitted intercept in log(dQ/dt) vs log(Q) plot')
    title('Change in recession intercept with season')
    fig_handles.SeasonalVarRecessions = fig;
end

end

