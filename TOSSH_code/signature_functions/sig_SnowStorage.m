function [SnowStorage, error_flag, error_str, fig_handles] = ...
    sig_SnowStorage(Q, t, P, varargin)
%sig_SnowStorage Calculates total snow storage using mass curve.
%   Difference in winter mass curves as a proxy for snow storage following
%   Schaefli (2016) and Horner et al. (2020).
%
%   INPUT
%   Q: streamflow [mm/timestep]
%   t: time [Matlab datetime]
%   P: precipitation [mm/timestep]
%   OPTIONAL
%   month_range: months over which the difference in the mass curves are
%       calculated, default = [10:12,1:3] (October to March)
%   plot_results: whether to plot results, default = false
%
%   OUTPUT
%   SnowStorage: estimated snow storage [mm]
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
%   P = data.P;
%   SnowStorage = sig_SnowStorage(Q,t,P);
%   SnowStorage = sig_SnowStorage(Q,t,P,'plot_results',true);
%
%   References
%   Schaefli, B., 2016. Snow hydrology signatures for model identification
%   within a limits of acceptability approach. Hydrological Processes,
%   30(22), pp.4019-4035.
%   Horner, I., Branger, F., McMillan, H., Vannier, O. and Braud, I.,
%   2020. Information content of snow hydrological signatures based on
%   streamflow, precipitation and air temperature. Hydrological Processes.
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
% time series have to be numeric and either a (n,1) or a (1,n) vector
addRequired(ip, 'Q', @(Q) isnumeric(Q) && (size(Q,1)==1 || size(Q,2)==1))
% date time series has to be numeric or datetime and either a (n,1) or a (1,n) vector
addRequired(ip, 't', @(t) (isnumeric(t) || isdatetime(t)) && (size(t,1)==1 || size(t,2)==1))
% time series have to be numeric and either a (n,1) or a (1,n) vector
addRequired(ip, 'P', @(P) isnumeric(P) && (size(P,1)==1 || size(P,2)==1))

% optional input arguments
% months has to be numeric and either a (n,1) or a (1,n) vector (default: winter months northern hemisphere)
addParameter(ip, 'month_range', [10:12,1:3], ...
    @(month_range) isnumeric(month_range) && (size(month_range,1)==1 || size(month_range,2)==1))
addParameter(ip, 'plot_results', false, @islogical) % whether to plot results

parse(ip, Q, t, P, varargin{:})
month_range = ip.Results.month_range;
plot_results = ip.Results.plot_results;

% create empty figure handle
fig_handles = [];

% data checks
[error_flag, error_str, timestep, t] = util_DataCheck(Q, t, 'P', P);
if error_flag == 2
    SnowStorage = NaN;
    return
end
timestep_days = days(timestep); % adjust for timestep

if any(month_range > 12) || any(month_range<1)
    error('month_range has to consist of values between 1 and 12.')
elseif ~all(diff(month_range) == 1 | diff(month_range) == -11)
    error('month_range has to consist of consecutive months.')
end

% calculate signature
% get average year
[Q_avg,t_avg] = util_AverageYear(Q,t,'start_water_year',month_range(1));
[P_avg,~] = util_AverageYear(P,t,'start_water_year',month_range(1));

% calculate maximum difference between cumulative sums
P_cumsum = cumsum(P_avg,'omitnan')./timestep_days;
Q_cumsum = cumsum(Q_avg,'omitnan')./timestep_days;
PQ_diff = P_cumsum - Q_cumsum;
is_month = ismember(month(t_avg),month_range);
[SnowStorage, index] = max(PQ_diff(is_month));

% TODO: other approach using inflection point in Horner et al. (2020)

% optional plotting
if plot_results
    day_of_year = 1:365;
    fig = figure('Position',[100 100 350 300]); hold on
    area(day_of_year(is_month),P_cumsum(diff(is_month)==-1)+0.*day_of_year(is_month),...
        'basevalue',0,'FaceColor',[.8 .8 .8],'FaceAlpha',0.2,'EdgeColor','none');
    plot(day_of_year,P_cumsum,'b','linewidth',1.5)
    plot(day_of_year,Q_cumsum,'r','linewidth',1.5)
    plot([index,index],[P_cumsum(index),Q_cumsum(index)],'k x-','linewidth',1.0)
    xlabel(strcat('Day of water year (starting in month',{' '},num2str(month_range(1)),')'))
    ylabel('Cumulative flow [mm]')
    legend('Considered months','P mass curve','Q mass curve','estimated snow storage','location','nw');
    fig_handles.SnowStorage = fig;
end

end
