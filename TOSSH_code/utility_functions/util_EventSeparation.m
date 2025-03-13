function [stormarray, error_flag, error_str, fig_handles] = util_EventSeparation(...
    dates, P, timestep, min_termination, min_duration, ...
    min_intensity_hour, min_intensity_day, ...
    min_intensity_hour_during, min_intensity_day_during, ...
    max_recessiondays, plot_results)
%util_EventSeparation takes rainfall data and picks out storm periods.
%
%   INPUT
%   dates: column vector of datenums
%   P: precipitation [mm/timestep]
%   timestep: time step of precipitation array [hours] (1=hourly, 24=daily)
%   min_termination: minimum termination time (time between storms) [hours]
%   min_duration: minimum duration of storm [hours]
%   min_intensity_hour: minimum intensity (per hour)
%   min_intensity_day: minimum intensity (per day)
%   min_intensity_hour_during: minimum timestep intensity allowed during
%       storm event without contributing to termination time
%   min_intensity_day_during: minimum timestep intensity allowed during
%       storm event without contributing to termination time
%   max_recessiondays: maximum number of days to allow recession after rain
%       ends
%   plot_results: whether to plot results, default = false
%
%   OUTPUT
%   stormarray: 2-column array with start and end locations of each storm
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
%   stormarray = util_EventSeparation(...
%       datenum(t), P, 1, 8, 5, 2, 10, 0.2, 1, 1, true);
%
%   Copyright (C) 2020
%   This software is distributed under the GNU Public License Version 3.
%   See <https://www.gnu.org/licenses/gpl-3.0.en.html> for details.

% check input parameters
if nargin < 11
    error('Not enough input arguments.')
end

% create empty figure handle
fig_handles = [];

% data checks

% default setting reads as good data
error_flag = 0;
error_str = '';

if ~ismember(timestep,[0.25,1,24])
    warning('Caution: The event separation function was designed for timesteps of 15 min, 1 hour or 1 day.')
end

if and(timestep<1,1/timestep ~= floor(1/timestep))
    error('Timestep must divide into 1 hour.')
end

% P is rainfall with timestep (hr) - 15 min (0.25), hour (1), day (24)

% create moving average series to check hourly/daily intensities
% create hourly moving average series
if timestep == 1
    P_hr = P;
elseif timestep < 1
    % calculate size of moving average window for 1 hr
    hr_window = 1/timestep;
    % append zeros, filter, remove zeros to center filter on each timestep
    P_hr = [P; zeros(floor(hr_window/2),1)];
    P_hr = filter((1/hr_window)*ones(1,hr_window),1,P_hr);
    P_hr = P_hr(1+floor(hr_window/2):end);
end
% create daily moving average series
if timestep == 24
    P_day = P;
elseif timestep < 24
    day_window = 24/timestep;
    P_day = [P; zeros(floor(day_window/2),1)];
    P_day = filter((1/day_window)*ones(1,day_window),1,P_day);
    P_day = P_day(1+floor(day_window/2):end);
end

% find gaps between storm events
% storm gaps when hourly rainfall below threshold for time greater than min_termination
if timestep <= 1
    % find all timesteps with hourly rainfall below threshold
    P_lowrain = P_hr <= min_intensity_hour_during;
else
    P_lowrain = P_day <= min_intensity_day_during;
end
P_lowrain(1) = 0;
% find beginning and end of runs of hourly rainfall below threshold
P_lowrain_change = P_lowrain(2:end)-P_lowrain(1:end-1);
begin_gap = find(P_lowrain_change == 1)+1;
end_gap = find(P_lowrain_change == -1);
% get complete gaps only
begin_gap = begin_gap(1:length(end_gap));
% get length of gaps in hours
length_gap = end_gap-begin_gap+1;
% identify too short gaps
short_gaps = find(length_gap < min_termination/timestep);
for i = 1:length(short_gaps) % delete these short gaps
    P_lowrain(begin_gap(short_gaps(i)):end_gap(short_gaps(i)))=0;
end

% check if potential storm periods meet criteria
% get potential storm periods
% all timesteps where intensity is high enough at hourly/daily timescale
potential_storms = (P_lowrain == 0);
potential_storms(1) = 0;
% identify runs (consecutive timesteps of rainfall intensity)
potential_storms_change = potential_storms(2:end)-potential_storms(1:end-1);
% identify beginning and end of storm periods
begin_storm = find(potential_storms_change == 1)+1;
end_storm = find(potential_storms_change == -1);
% remove last 'beginning of storm' if it does not complete within time series
begin_storm = begin_storm(1:length(end_storm));

valid_storm = zeros(length(begin_storm),3);
% cycle through potential storms and check if valid
for i = 1:length(begin_storm)
    
    % check duration
    if end_storm(i) - begin_storm(i) + 1 >= min_duration/timestep
        valid_storm(i,1) = 1;
    end
    
    % check hourly intensity (if timestep <= 1 hr)
    if timestep <= 1
        if max(P_hr(begin_storm(i):end_storm(i))) >= min_intensity_hour
            valid_storm(i,2) = 1;
        end
    end
    
    % check daily intensity
    if max(P_day(begin_storm(i):end_storm(i))) >= min_intensity_day/(24/timestep)
        valid_storm(i,3) = 1;
    end
end

% valid storm should have long enough duration and high enough intensity at either hourly or daily timescale
valid_overall = and(valid_storm(:,1),or(valid_storm(:,2),valid_storm(:,3)));

% oputput array records beginning and end of storm periods
stormarray = [begin_storm(valid_overall),end_storm(valid_overall)];
% disp(['Number of storm events: ', num2str(size(stormarray,1))]);

if size(stormarray,1)==0
    error_flag = 3;
    error_str = ['Error: No events detected. ', error_str];
    return
elseif size(stormarray,1) < 10
    error_flag = 1;
    error_str = ['Warning: Fewer than 10 events detected, results might not be robust. ', error_str];
end

% Get suitable end of storm event for signatures that use flow. An event 
% goes for either 5 days after rain, or until the next event starts, or
% until rainfall is greater than the min_intensity_hour_during.

% max-day criterion or until next event
stormarray(1:end-1,3) = min(stormarray(1:end-1,2)+max_recessiondays*24/timestep,stormarray(2:end,2)-1);
stormarray(end,3) = min(stormarray(end,2)+max_recessiondays*24/timestep,length(P));

% until rainfall is over the maximum
if timestep <= 1
    % find all timesteps with hourly rainfall below threshold
    P_lowrain = P_hr > min_intensity_hour_during;
else
    P_lowrain = P_day > min_intensity_day_during;
end

% find length of recession period after rainfall ends, before next storm begins
recession_rain = zeros(size(stormarray,1),1);
for i = 1:size(stormarray,1)
    rain_index = find(P_lowrain((stormarray(i,2)+1):stormarray(i,3))==1,1,'first');
    if numel(rain_index) > 0
        recession_rain(i) = rain_index-1;
    else
        recession_rain(i) = inf;
    end
end
stormarray(:,3) = min(stormarray(:,3),stormarray(:,2)+recession_rain);

% optional plotting
if plot_results
    % if plotting requested, show rainfall with overlaid storm events
    dates_dt = datetime(dates,'ConvertFrom','datenum');
    fig = figure('Position',[100 100 700 250]);
    P_max = max(P);
    hold on
    for i = 1:size(stormarray,1)
        h1=fill([dates_dt(stormarray(i,1)),dates_dt(stormarray(i,1)),...
            dates_dt(stormarray(i,2)),dates_dt(stormarray(i,2))],...
            [0, P_max, P_max, 0],[106, 174, 214]/255,'LineStyle','none');
        h2=fill([dates_dt(stormarray(i,2)),dates_dt(stormarray(i,2)),...
            dates_dt(stormarray(i,3)),dates_dt(stormarray(i,3))],...
            [0, P_max, P_max, 0],[198, 219, 239]/255,'LineStyle','none');
    end
    h3=plot(dates_dt(:),P(:),'k-','linewidth',1.0);
    xlabel('Date')
    ylabel('Rainfall [mm/timestep]')
    legend([h3,h1,h2],{'Rainfall','Storm Period','Recession Period'})
    fig_handles.EventSeparation = fig;
    
end
