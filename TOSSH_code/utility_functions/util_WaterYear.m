function [water_year] = util_WaterYear(t, varargin)
%util_WaterYear gets the water year associated with a date.
%
%   INPUT
%   t: time [Matlab datetime, scalar or array]
%   OPTIONAL
%   WY_start_month: numeric month in which water year starts
%
%   OUTPUT
%   water_year: water year in which date falls
%
%   EXAMPLE
%   % load example data 
%   data = load('example/example_data/33029_daily.mat'); 
%   t = data.t;
%   water_year = util_WaterYear(t, 'WY_start_month',10);
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

% optional input arguments
addParameter(ip, 'WY_start_month', 10, ...
    @(WY_start_month) (isnumeric(WY_start_month)&&floor(WY_start_month)==WY_start_month) ...
    && (WY_start_month>=1 && WY_start_month<=12)) 

parse(ip, t, varargin{:})
WY_start_month = ip.Results.WY_start_month;

% % timestep checks
% if isnumeric(t)
%     t = datetime(t,'ConvertFrom','datenum');
%     warning('Converted datenum to datetime.')
% end

% get month and year of each datetime
water_year = year(t);
month_vals = month(t);

% get months associated with previous year and subtract 1 from year value
% water_year(month_vals<WY_start_month) = water_year(month_vals<WY_start_month) - 1;  

% get months associated with previous year and add 1 from year value 
% USGS: the water year is designated by the calendar year in which it ends
water_year(month_vals>WY_start_month) = water_year(month_vals>WY_start_month) + 1;     

end
