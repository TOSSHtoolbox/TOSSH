function [EventRR, error_flag, error_str, fig_handles] = sig_EventRR(Q, t, P, varargin)
%sig_EventRR calculates event runoff ratio.
%   Extracts events from the rainfall time series and calculates the
%   average runoff ratio over all events.
%
%   INPUT
%   Q: streamflow [mm/timestep]
%   t: time [Matlab datetime]
%   P: precipitation [mm/timestep]
%
%   OPTIONAL
%   (All these parameters enable the user to tweak the definition of storm
%   events:)
%   min_termination: minimum termination time between storm events [hours],
%       default = 8
%   min_duration: minimum duration of storm [hours], default = 5
%   min_intensity_hour: minimum hourly rainfall [mm/hour], default = 2
%   min_intensity_day: minimum daily rainfall [mm/day], default = 10
%   min_intensity_hour_during: minimum rainfall allowed during an event
%       without contributing to termination time [mm/hour], default = 0.2
%   min_intensity_day_during: minimum rainfall allowed during an event
%       without contributing to termination time [mm/day], default = 1
%   max_recessiondays: maximum length of recession period after rainfall
%       [days] deemed to be part of the preceeding event (shorter if
%       another event starts), default = 1
%   plot_results: whether to plot results, default = false
%
%   OUTPUT
%   EventRR: event runoff ratio [-]
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
%   EventRR = sig_EventRR(Q, t, P);
%   EventRR = sig_EventRR(Q, t, P, 'plot_results',true);
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
% flow time series have to be numeric and either a (n,1) or a (1,n) vector
addRequired(ip, 'Q', @(Q) isnumeric(Q) && (size(Q,1)==1 || size(Q,2)==1))
% date time series has to be numeric or datetime and either a (n,1) or a (1,n) vector
addRequired(ip, 't', @(t) ((isnumeric(t) || isdatetime(t))) && ((size(t,1)==1 || size(t,2)==1)))
% precipitation time series have to be numeric and either a (n,1) or a (1,n) vector
addRequired(ip, 'P', @(P) isnumeric(P) && (size(P,1)==1 || size(P,2)==1))

% optional input arguments
addParameter(ip, 'min_termination', 8, @isnumeric) % minimum termination time
% (time between storms) in hours
addParameter(ip, 'min_duration', 5, @isnumeric) % minimum duration of storm in hour
addParameter(ip, 'min_intensity_hour', 2, @isnumeric) % minimum hourly rainfall (per hour)
addParameter(ip, 'min_intensity_day', 10, @isnumeric) % minimum daily rainfall (per day)
addParameter(ip, 'min_intensity_hour_during', 0.2, @isnumeric) % minimum timestep
% intensity allowed during storm event without contributing to termination time
addParameter(ip, 'min_intensity_day_during', 1, @isnumeric) % minimum timestep
% intensity allowed during storm event without contributing to termination time
addParameter(ip, 'max_recessiondays', 1, @isnumeric) % maximum number of days to allow recession after rain
addParameter(ip, 'plot_results', false, @islogical) % whether to plot results (1 graph)

parse(ip, Q, t, P, varargin{:})

min_termination = ip.Results.min_termination;
min_duration = ip.Results.min_duration;
min_intensity_hour = ip.Results.min_intensity_hour;
min_intensity_day = ip.Results.min_intensity_day;
min_intensity_hour_during = ip.Results.min_intensity_hour_during;
min_intensity_day_during = ip.Results.min_intensity_day_during;
max_recessiondays = ip.Results.max_recessiondays;
plot_results = ip.Results.plot_results;

% create empty figure handle
fig_handles = [];

% data checks
[error_flag, error_str, timestep, t] = util_DataCheck(Q, t, 'P', P);
if error_flag == 2
    EventRR = NaN;
    return
end
timestep = hours(timestep);

% run event separation algorithm
error_flag_tmp = error_flag; % temporarily store error flag from data check
error_str_tmp = error_str;
[stormarray, error_flag, error_str, fig_handles] = util_EventSeparation(...
    datenum(t), P, timestep, min_termination, min_duration, ...
    min_intensity_hour, min_intensity_day, ...
    min_intensity_hour_during, min_intensity_day_during, ...
    max_recessiondays, plot_results);
if error_flag == 3
    EventRR = NaN;
    return
else
    error_flag = max([error_flag_tmp, error_flag]);
    error_str = [error_str_tmp, error_str];
end

% cycle through storms and sum rainfall and flow for each one, calculate
% runoff ratio
event_rr_array = zeros(size(stormarray));
for i = 1:size(stormarray,1)
    event_rr_array(i,1) = sum(P(stormarray(i,1):stormarray(i,3)));
    event_rr_array(i,2) = sum(Q(stormarray(i,1):stormarray(i,3)));
    event_rr_array(i,3) = event_rr_array(i,2)/event_rr_array(i,1);
end

% return average event runoff ratio
EventRR = mean(event_rr_array(:,3),'omitnan');

end