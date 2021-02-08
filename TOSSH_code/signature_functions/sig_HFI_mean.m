function [HFI_mean, error_flag, error_str] = sig_HFI_mean(Q,t, varargin)
%sig_HFI_mean calculates mean half flow interval.
%   Calculates time span between the date on which the cumulative discharge 
%   since start of water year (default: October) reaches (here: exceeds) a 
%   quarter of the annual discharge and the date on which the cumulative 
%   discharge since start of water year (default: October) reaches three 
%   quarters of the annual discharge.
%
%   INPUT
%   Q: streamflow [mm/timestep]
%   t: time [Matlab datenum]
%   OPTIONAL
%   start_water_year: first month of water year, default = 10 (October)
%
%   OUTPUT
%   HFI_mean: mean half flow interval [days]
%   error_flag: 0 (no error), 1 (warning), 2 (error in data check), 3
%       (error in signature calculation)
%   error_str: string contraining error description
%
%   EXAMPLE
%   % load example data 
%   data = load('example/example_data/33029_daily.mat'); 
%   Q = data.Q; 
%   t = data.t;
%   HFI_mean = sig_HFI_mean(Q,t);
%   HFI_mean = sig_HFI_mean(Q,t,'start_water_year',1);
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
validationFcn = @(x) isnumeric(x) && isscalar(x) && (x >= 1) && (x <= 12) && floor(x)==x;
addParameter(ip, 'start_water_year', 10, validationFcn) % when does the water year start? Default: 10

parse(ip, Q, t, varargin{:})
start_water_year = ip.Results.start_water_year;

% data checks
[error_flag, error_str, timestep, t] = util_DataCheck(Q, t);
if error_flag == 2
    HFI_mean = NaN;
    return
end
timestep_days = days(timestep); % adjust for timestep

% calculate signature
% get individual years
[year_vec, month_vec, day_vec] = ymd(t);
year_start = min(year_vec);
year_end = max(year_vec);
year_list = [year_start:year_end]';

Q_temp = Q;
% Q_annual = NaN(year_end-year_start,1);
% Q_daily = NaN(365,year_end-year_start);
HFI = NaN(year_end-year_start,1);

% extract years
error_tmp = false;
for y = 2:length(year_list) % since we use water years, we always start in the "2nd year"
    year = year_list(y);
    Q_water_year = ...
        [Q_temp(year_vec==year-1 & month_vec>=start_water_year); ...
        Q_temp(year_vec==year & month_vec<start_water_year)];
    Q_25_sum = 0.25*sum(Q_water_year);
    Q_75_sum = 0.75*sum(Q_water_year);
    Q_cumsum = cumsum(Q_water_year);
    aux_index = 1:length(Q_water_year);
    HFI_aux_25 = aux_index(Q_cumsum>Q_25_sum);
    HFI_aux_75 = aux_index(Q_cumsum>Q_75_sum);    
    if isempty(HFI_aux_25) || isempty(HFI_aux_75) % if there is no flow          
        error_tmp = true;        
    else
        HFI_25 = HFI_aux_25(1);
        HFI_75 = HFI_aux_75(1);
        HFI(y-1) = HFI_75 - HFI_25;
    end
end

if error_tmp
    error_flag = 1;
    error_str = ['Warning: Some years have no flow. ', error_str];
end

% get mean half flow interval (ignoring NaNs)
HFI_mean = mean(HFI,'omitnan')*timestep_days;

end
