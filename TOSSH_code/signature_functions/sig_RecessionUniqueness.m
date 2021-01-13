function [Spearmans_rho, error_flag, error_str, fig_handles] ...
    = sig_RecessionUniqueness(Q, t, varargin)
%sig_RecessionUniqueness calculates uniqueness of storage-discharge relationship.
%   Calculates Spearman's rank correlation between Q and dQ/dt to measure 
%   how unique the storage-discharge relationship is.
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
%   dQdt_method: method for dQ/dt calculation, default = 'ETS'
%   plot_results: whether to plot results, default = false
%
%   OUTPUT
%   Spearmans_rho: Spearman's rank correlation between Q and dQ/dt
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
%   Spearmans_rho = sig_RecessionUniqueness(Q,t,'plot_results',true);
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
addParameter(ip, 'ignoreNaN', 'y', @ischar) % ignore NaN values y/n?
addParameter(ip, 'recession_length', 5, @isnumeric) % length of decreasing 
% flow section (amount of timesteps) to be declared a recession
addParameter(ip, 'n_start', 1, @isnumeric) % days to be removed at beginning of recession
addParameter(ip, 'eps', 0, @isnumeric) % allowed increase in flow during recession period
addParameter(ip, 'start_of_recession', 'peak', @ischar) % defines start of a recession
addParameter(ip, 'filter_par', 0.925, @isnumeric) % smoothing parameter of 
% Lyne-Hollick Filter to determine start of recession (higher = later recession start)
addParameter(ip, 'dQdt_method', 'ETS', @ischar) % how to calculate dQ/dt
addParameter(ip, 'plot_results', false, @islogical) % whether to plot results (2 graphs)

parse(ip, Q, t, varargin{:})
recession_length = ip.Results.recession_length;
n_start = ip.Results.n_start;
eps = ip.Results.eps;
start_of_recession = ip.Results.start_of_recession;
filter_par = ip.Results.filter_par;
dQdt_method = ip.Results.dQdt_method;
plot_results = ip.Results.plot_results;

% create empty figure handle
fig_handles = [];

% data checks
[error_flag, error_str, timestep, t] = util_DataCheck(Q, t);
if error_flag == 2
    Spearmans_rho = NaN;
    return
end

% calculate signature

% get recession segments
error_flag_tmp = error_flag; % temporarily store error flag from data check
error_str_tmp = error_str;
[flow_section, error_flag, error_str, fig_handles] = util_RecessionSegments(Q, t, ...
    'recession_length', recession_length, 'eps', eps, ...
    'filter_par', filter_par, 'plot_results', plot_results, ...
    'start_of_recession', start_of_recession, 'n_start', n_start);
if error_flag == 3
    Spearmans_rho = NaN;
    return
else
    error_flag = max([error_flag_tmp, error_flag]);
    error_str = [error_str_tmp, error_str];
end

% get flow rate gradient and corresponding flows
[dQdt, Qm, ~, ~] = ...
    util_dQdt(Q, t, flow_section, 'method', dQdt_method);

% calculate Spearman's rho between Q and dQ/dt
rec = ~isnan(Qm);
Spearmans_rho = corr(Qm(rec),dQdt(rec),'Type','Spearman');

if plot_results
    fig = figure('Position',[100 100 350 300]); hold on
    plot(Qm,-dQdt,'k .','linewidth',2)
    xlabel('Q [mm/timestep]') 
    ylabel('-dQ/dt [mm/timestep^2]') 
    set(gca,'XScale','log')
    set(gca,'YScale','log')
    title(['Spearmans rho:','{ }',num2str(Spearmans_rho)])
    fig_handles.RecessionUniqueness = fig;
end

end

