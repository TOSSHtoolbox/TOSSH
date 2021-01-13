function [X_avg, t_avg, fig_handles] = util_AverageYear(X, t, varargin)
%util_AverageYear calculates average year, i.e. the mean on each day.
%   Note that the mean is returned in mm/timestep, i.e. it has to be 
%   converted to mm/day if the daily flow is required. 
%
%   INPUT
%   X: time series, e.g. streamflow [mm/timestep]
%   t: time [Matlab datetime]
%   OPTIONAL
%   start_water_year: first month of water year, default = 1 (January)
%   plot_results: whether to plot results, default = false
%
%   OUTPUT
%   X_avg: average flow on each date [mm/timestep]
%   t_avg: corresponding dates
%   fig_handles: figure handles to manipulate figures (empty if plotting is
%       not requested)
%
%   EXAMPLE
%   % load example data
%   data = load('example/example_data/33029_daily.mat');
%   Q = data.Q;
%   t = data.t;
%   Q_avg = util_AverageYear(Q,t);
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
addRequired(ip, 'X', @(X) isnumeric(X) && (size(X,1)==1 || size(X,2)==1))
% date time series has to be numeric or datetime and either a (n,1) or a (1,n) vector
addRequired(ip, 't', @(t) (isnumeric(t) || isdatetime(t)) && (size(t,1)==1 || size(t,2)==1))

% optional input arguments
addParameter(ip, 'start_water_year', 1, @isnumeric)
addParameter(ip, 'plot_results', false, @islogical) % whether to plot results 

parse(ip, X, t, varargin{:})
start_water_year = ip.Results.start_water_year;
plot_results = ip.Results.plot_results;

% create empty figure handle
fig_handles = [];

% calculate average year
months = [start_water_year:12, 1:start_water_year-1]';
n_days = [31 28 31 30 31 30 31 31 30 31 30 31]';

X_avg = NaN(365,1);
t_avg = NaN(365,3);
index = 1;

% loop over months
day_vec = day(t);
month_vec = month(t);
for i = 1:12
    m = months(i);
    % loop over days
    for d = 1:n_days(m)
        % X_avg(index) = mean(X(day_vec == d & month_vec == m),'omitnan');
        X_tmp = X(day_vec == d & month_vec == m);
        X_avg(index) = sum(X_tmp,'omitnan')/length(X_tmp);
        t_avg(index,:) = [0,m,d]; % dummy year
        index = index + 1;
    end
end
t_avg = datetime(t_avg); % dummy year

% optional plotting
if plot_results
    fig = figure('pos',[100 100 350 300]); hold on
    plot(t_avg,X_avg,'.')
    xlabel('Day')
    ylabel('Mean on that day [mm/timestep]')
    fig_handles.AverageYear = fig;
end

end
