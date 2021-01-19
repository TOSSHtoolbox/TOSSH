function [X_annual, X_monthly, year_list, days_per_year, error_flag, error_str] = ...
    util_AggregateTimeSeries(X, t, start_water_year)
%util_AggregateTimeSeries calculates annual and monthly sums of time series.
%   Start of water year can be specified so that the average is taken over
%   a water year.
%
%   INPUT
%   X: time series - should be continuous (e.g. Q, typically [mm/timestep])
%   t: time [Matlab datetime]
%   start_water_year: first month of water year, default = 1 (January)
%
%   OUTPUT
%   X_annual: annual sums [mm/year]
%   X_monthly: monthly sums [mm/month]
%   year_list: corresponding years
%   days_per_year: amount of non-NaN days per year
%   error_flag: 0 (no error), 1 (warning), 2 (error in data check), 3
%       (error in signature calculation)
%   error_str: string contraining error description
%
%   EXAMPLE
%   % load example data
%   data = load('example/example_data/33029_daily.mat');
%   Q = data.Q;
%   t = data.t;
%   start_water_year = 10;
%   [Q_annual, Q_monthly, year_list, leap_year] = ...
%       util_AggregateTimeSeries(Q, t, start_water_year);
%
%   Copyright (C) 2020
%   This software is distributed under the GNU Public License Version 3.
%   See <https://www.gnu.org/licenses/gpl-3.0.en.html> for details.

if nargin < 3
    start_water_year = 1;
end

% default setting reads as good data
error_flag = 0;
error_str = '';

% get years and months
[year_vec, month_vec, day_vec] = ymd(t);

if start_water_year == 1
    % calendar year
    year_start = min(year_vec);
    year_end = max(year_vec);
else
    % water year always corresponds to the last day of the water year,
    % e.g. water year starting from 1 October 1999 is water year 2000
    year_start = min(year_vec)+1;
    year_end = max(year_vec);
end
year_list = [year_start:year_end]';

if month_vec(1) ~= start_water_year && day_vec(1)~=1
    error_flag = 1;
    error_str = ['Warning: Time series and water year do not match. Incomplete years possible. ', error_str];
elseif month_vec(end) ~= start_water_year-1 && day_vec(end)<28
    error_flag = 1;
    error_str = ['Warning: Time series and water year do not match. Incomplete years possible. ', error_str];
end

X_annual = NaN(year_end-year_start,1);
X_monthly = NaN(year_end-year_start,12);
days_per_year = NaN(size(year_list));

% extract years and months
for y = 1:length(year_list)
    year = year_list(y);
    X_water_year = ...
        [X(year_vec==year-1 & month_vec>=start_water_year); ...
        X(year_vec==year & month_vec<start_water_year)];
    X_annual(y) = sum(X_water_year,'omitnan');
    days_per_year(y) = length(X_water_year(~isnan(X_water_year)));        
    for m = start_water_year:12
        X_monthly(y,m-start_water_year+1) = sum(X(year_vec==year-1 & month_vec==m),'omitnan');
    end
    for m = 1:start_water_year-1
        X_monthly(y,m+13-start_water_year) = sum(X(year_vec==year & month_vec==m),'omitnan');
    end
end

end