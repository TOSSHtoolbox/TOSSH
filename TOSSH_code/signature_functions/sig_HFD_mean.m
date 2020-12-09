function [HFD_mean, error_flag, error_str] = sig_HFD_mean(Q, t, varargin)
%sig_HFD_mean calculates mean half flow date.
%   Calculates day since start of water year on which the cumulative 
%   discharge (default: October) reaches half of the annual discharge.
%
%   INPUT
%   Q: streamflow [mm/timestep]
%   t: time [Matlab datetime]
%   OPTIONAL
%   start_month: starting month, default = 10 (October)
%
%   OUTPUT
%   HFD_mean: mean half flow date [day since start of water year]
%   error_flag: 0 (no error), 1 (warning), 2 (error in data check), 3
%       (error in signature calculation)
%   error_str: string contraining error description
%
%   EXAMPLE
%   % load example data 
%   data = load('example/example_data/33029_daily.mat'); 
%   Q = data.Q; 
%   t = data.t;
%   HFD_mean = sig_HFD_mean(Q,t);
%   HFD_mean = sig_HFD_mean(Q,t,'start_month',1);
%
%   References
%   Court, A., 1962. Measures of streamflow timing. Journal of Geophysical
%   Research, 67(11), pp.4335-4339.
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
addParameter(ip, 'start_month', 10, @(month) isnumeric(month) && numel(month)==1) 

parse(ip, Q, t, varargin{:})
start_month = ip.Results.start_month;

% data checks
[error_flag, error_str, timestep, t] = util_DataCheck(Q, t);
if error_flag == 2
    HFD_mean = NaN;
    return
end
timestep_days = days(timestep); % adjust for timestep

if any(start_month<1 | start_month>12) || any(floor(start_month)~=start_month)
    error('Month has to be a vector containing integers between 1 and 12.')
end

% calculate signature
% get individual years
[year_vec, month_vec, day_vec] = ymd(t);
year_start = min(year_vec);
year_end = max(year_vec);
year_list = [year_start:year_end]';

Q_temp = Q;
% Q_annual = NaN(year_end-year_start,1);
% Q_daily = NaN(365,year_end-year_start);
HFD = NaN(year_end-year_start,1);

% extract years
error_tmp = false;
for y = 2:length(year_list) % since we use water years, we always start in the "2nd year"
    try
        year = year_list(y);
        Q_water_year = ...
            [Q_temp(year_vec==year-1 & month_vec>=start_month); ...
            Q_temp(year_vec==year & month_vec<start_month)];
        Q_half_sum = 0.5*sum(Q_water_year);
        Q_cumsum = cumsum(Q_water_year);
        aux_index = 1:length(Q_water_year);
        HFD_aux = aux_index(Q_cumsum>=Q_half_sum);
        HFD(y-1) = HFD_aux(1);
    catch
        error_tmp = true;
    end
end

if error_tmp
    error_flag = 1;
    error_str = ['Warning: Years containing NaN values are ignored. ', error_str];
end

% get mean half flow date
HFD_mean = mean(HFD,'omitnan')*timestep_days;

end


