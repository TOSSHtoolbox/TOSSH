function [BaseflowRecessionK, error_flag, error_str, fig_handles] = ...
    sig_BaseflowRecessionK(Q, t, varargin)
%sig_BaseflowRecessionK calculates baseflow recession constant.
%   Calculates baseflow recession constant assuming exponential recession
%   behaviour (Safeeq et al., 2013). Master recession curve (MRC) is
%   constructed using the adapted matching strip method (Posavec et al.,
%   2006).
%
%   Notes:
%   According to Safeeq et al. (2013), K<0.065 represent groundwater
%   dominated slow-draining systems, K>=0.065 represent shallow subsurface
%   flow dominated fast draining systems.
%   (to do: remove snow-affected sections of the time series)
%
%   INPUT
%   Q: streamflow [mm/timestep]
%   t: time [Matlab datetime]
%   OPTIONAL
%   recession_length: min. length of recessions [days], default = 15
%   n_start: days to be removed after start of recession
%   eps: allowed increase in flow during recession period, default = 0
%   start_of_recession: define start of recession when baseflow filter
%       rejoins the curve ("baseflow"), or after hydrograph peak ("peak")
%   fit_method: method to fit MRC, default = 'nonparametric_analytic'
%   filter_par: smoothing parameter of Lyne-Hollick filter to determine
%      start of recession (higher = later recession start), default = 0.925
%   plot_results: whether to plot results, default = false
%
%   OUTPUT
%   BaseflowRecessionK: Baseflow recession constant [1/timestep]
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
%   BaseflowRecessionK = sig_BaseflowRecessionK(Q,t);
%   BaseflowRecessionK = sig_BaseflowRecessionK(Q,t,'plot_results',true,'recession_length',5);
%
%   References
%   Safeeq, M., Grant, G.E., Lewis, S.L. and Tague, C.L., 2013. Coupling
%   snowpack and groundwater dynamics to interpret historical streamflow
%   trends in the western United States. Hydrological Processes, 27(5),
%   pp.655-668.
%   Posavec, K., Bacani, A. and Nakic, Z., 2006. A visual basic spreadsheet
%   macro for recession curve analysis. Groundwater, 44(5), pp.764-767.
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
addParameter(ip, 'recession_length', 15, @isnumeric) % length of decreasing
% flow in days to be declared a recession
addParameter(ip, 'n_start', 0, @isnumeric) % days to be removed at beginning of recession
addParameter(ip, 'eps', 0, @isnumeric) % allowed increase in flow during recession period
addParameter(ip, 'start_of_recession', 'baseflow', @ischar) % defines start of a recession
addParameter(ip, 'fit_method', 'nonparametric_analytic', @ischar) % how to fit MRC
addParameter(ip, 'filter_par', 0.925, @isnumeric) % smoothing parameter of
% Lyne Hollick Filter to determine start of recession (higher = later recession start)
addParameter(ip, 'plot_results', false, @islogical) % whether to plot results (2 graphs)

parse(ip, Q, t, varargin{:})
recession_length = ip.Results.recession_length;
n_start = ip.Results.n_start;
eps = ip.Results.eps;
start_of_recession = ip.Results.start_of_recession;
fit_method = ip.Results.fit_method;
filter_par = ip.Results.filter_par;
plot_results = ip.Results.plot_results;

% create empty figure handle
fig_handles = [];

% data checks
[error_flag ,error_str, timestep, t] = util_DataCheck(Q, t);
if error_flag == 2
    BaseflowRecessionK = NaN;
    return
end

% calculate signature
% steps from Safeeq et al. (2013)

% identify all individual baseflow recession segments
error_flag_tmp = error_flag; % temporarily store error flag from data check
error_str_tmp = error_str;
[flow_section, error_flag, error_str, fig_handles] = util_RecessionSegments(Q, t, ...
    'recession_length', recession_length, 'eps', eps, ...
    'filter_par', filter_par, 'plot_results', plot_results, ...
    'start_of_recession', start_of_recession, 'n_start', n_start);
if error_flag == 3
    BaseflowRecessionK = NaN;
    return
else
    error_flag = max([error_flag_tmp, error_flag]);
    error_str = [error_str_tmp, error_str];
end

% From Safeeq et al. (2013): To minimize the effect of snowmelt recharge on
% k, recession segments identified between the onset of the
% snowmelt-derived streamflow pulse and 15 August were excluded. Days of
% snowmelt pulse onset were determined following the method of Cayan et al.
% (2001). This is not implemented. Requires additional information on the
% time series such as whether it is in a snow region, and requires series
% to consist of complete years only.

% MRC constucted using the adapted matching strip method (Posavec et al., 2006)
[MRC] = util_MasterRecessionCurve(Q, flow_section,'fit_method',fit_method,'match_method','log','plot_results',false);

% k = slope of the linear regression between log-transformed discharge and recession length
mdl = [(MRC(:,1).^0) (MRC(:,1))]\log(MRC(:,2));
BaseflowRecessionK = -mdl(2);
if ~isreal(BaseflowRecessionK)
    error_flag = 3;
    error_str = ['Error: Complex BaseflowRecessionK. ', error_str];
    BaseflowRecessionK = NaN;
    return
end

%If timestep is not daily, convert units of BaseflowRecessionK to 1/day
timestep_factor = days(1)/timestep;
BaseflowRecessionK = BaseflowRecessionK * timestep_factor;

% optional plotting
if plot_results
    fig = figure('Position',[100 100 750 300]);
    subplot(1,2,1)
    plot(MRC(:,1),(MRC(:,2)),'kx')
    hold on
    plot(sort(MRC(:,1)),exp(mdl(1) + mdl(2).*(sort(MRC(:,1)))),'g-','linewidth',2)
    xlabel('Relative time [timestep]')
    ylabel('Flow [mm/timestep]')
    title('Fitted recession curve')
    legend('MRC','Exponential fit')
    
    subplot(1,2,2)
    semilogy(MRC(:,1),(MRC(:,2)),'kx')
    hold on
    semilogy(sort(MRC(:,1)),exp(mdl(1) + mdl(2).*(sort(MRC(:,1)))),'g-','linewidth',2)
    xlabel('Relative time [timestep]')
    ylabel('Flow [mm/timestep]')
    title('Fitted recession curve - log scale')
    legend('MRC','Exponential fit')
    
    fig_handles.BaseflowRecessionK = fig;
end

end

