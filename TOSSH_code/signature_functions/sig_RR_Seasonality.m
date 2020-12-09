function [RR_Seasonality, error_flag, error_str] = sig_RR_Seasonality(Q, t, P, varargin)
%sig_RR_Seasonality calculates seasonality of runoff ratio.
%   Ratio between summer and winter runoff ratio (i.e. fraction of 
%   precipitation that leaves the watershed as flow).
%
%   INPUT
%   Q: streamflow [mm/timestep]
%   t: time [Matlab datetime]
%   P: precipitation [mm/timestep]
%   OPTIONAL
%   summer_start: month when 6-month summer is deemed to start, default = 4
%       (April)
%
%   OUTPUT
%   RR_Seasonality: ratio between summer and winter runoff ratios [-]
%   error_flag: 0 (no error), 1 (warning), 2 (error in data check), 3
%       (error in signature calculation)
%   error_str: string contraining error description
%
%   EXAMPLE
%   % load example data 
%   data = load('example/example_data/33029_daily.mat'); 
%   Q = data.Q; 
%   t = data.t;
%   P = data.P;
%   % example for northern hemisphere when summer season starts in April
%   RR_seasonality = sig_RR_Seasonality(Q, t, P, 'summer_start', 4);
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
% time series have to be numeric and either a (n,1) or a (1,n) vector
addRequired(ip, 'P', @(P) isnumeric(P) && (size(P,1)==1 || size(P,2)==1)) 

% optional input arguments
addParameter(ip, 'summer_start', 4, @isnumeric) % month when 6-month summer is deemed to start

parse(ip, Q, t, P, varargin{:})
summer_start = ip.Results.summer_start;

% data checks
[error_flag, error_str, timestep, t] = util_DataCheck(Q, t, 'P', P);
if error_flag == 2
    RR_Seasonality = NaN;
    return
end

% calculate signature
% get summer timesteps
summer_months = mod([summer_start:summer_start+5],12);
summer_months(summer_months==0)=12;
summer_index = find(ismember(month(t),summer_months));
winter_index = find(~ismember(month(t),summer_months));

% get summer and winter runoff ratios
SummerRR = mean(Q(summer_index),'omitnan')./mean(P(summer_index),'omitnan');
WinterRR = mean(Q(winter_index),'omitnan')./mean(P(winter_index),'omitnan');

% ratio of summer and winter RRs
RR_Seasonality = SummerRR/WinterRR;

end
