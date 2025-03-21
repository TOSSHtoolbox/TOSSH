function [water_year] = util_WaterYear(t, varargin)
%util_WaterYear gets the water year associated with a date.
%
%   INPUT
%   t: time [Matlab datetime, scalar or array]
%   OPTIONAL
%   start_water_year: numeric month in which water year starts
%
%   OUTPUT
%   water_year: water year in which date falls
%
%   EXAMPLE
%   % load example data 
%   data = load('example/example_data/33029_daily.mat'); 
%   t = data.t;
%   water_year = util_WaterYear(t, 'start_water_year',10);
%
%   Copyright (C) 2020
%   This software is distributed under the GNU Public License Version 3.
%   See <https://www.gnu.org/licenses/gpl-3.0.en.html> for details.

% check input parameters
if nargin < 1
    error('Not enough input arguments.')
end

ip = inputParser;

% required input arguments
% date time series has to be numeric or datetime and either a (n,1) or a (1,n) vector
addRequired(ip, 't', @(t) (isnumeric(t) || isdatetime(t)) && (size(t,1)==1 || size(t,2)==1))

% optional input argument
validationFcn = @(x) isnumeric(x) && isscalar(x) && (x >= 1) && (x <= 12) && floor(x)==x;
addParameter(ip, 'start_water_year', 10, validationFcn) % when does the water year start?

parse(ip, t, varargin{:})
start_water_year = ip.Results.start_water_year;

% % timestep checks
% if isnumeric(t)
%     t = datetime(t,'ConvertFrom','datenum');
%     warning('Converted datenum to datetime.')
% end

% get month and year of each datetime
water_year = year(t);
month_vals = month(t);

% get months associated with previous year and subtract 1 from year value
% water_year(month_vals<start_water_year) = water_year(month_vals<start_water_year) - 1;  

% get months associated with previous year and add 1 from year value 
% USGS: the water year is designated by the calendar year in which it ends
water_year(month_vals>=start_water_year) = water_year(month_vals>=start_water_year) + 1;     

end
