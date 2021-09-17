function [error_flag, error_str, timestep, t] = util_DataCheck(Q, t, varargin)
%util_DataCheck checks data for various things.
%   Checks data for unrealistic values, NaN in time series, inconsistent
%   time series lenghts, etc., and returns warnings or errors.
%
%   INPUT
%   Q: streamflow [mm/timestep]
%   t: time [Matlab datetime]
%   OPTIONAL
%   P: precipitation [mm/timestep]
%   PET: potential evapotranspiration [mm/timestep]
%   T: temperature [degC]
%
%   OUTPUT
%   error_flag: 0 (no error), 1 (warning), 2 (error in data check), 3
%       (error in signature calculation)
%   error_str: string contraining error description
%   timestep: (median) timestep of data
%   t: time in Matlab datetime format
%
%   EXAMPLE
%   % load example data
%   data = load('example/example_data/33029_daily.mat');
%   Q = data.Q;
%   t = data.t;
%   P = data.P;
%   dataCheck = util_DataCheck(Q, t);
%   dataCheck = util_DataCheck(Q, t, 'P', P);
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
% P has to be numeric and either a (n,1) or a (1,n) vector
addParameter(ip, 'P', [], @(P) isnumeric(P) && (size(P,1)==1 || size(P,2)==1))
% PET has to be numeric and either a (n,1) or a (1,n) vector
addParameter(ip, 'PET', [], @(PET) isnumeric(PET) && (size(PET,1)==1 || size(PET,2)==1))
% T has to be numeric and either a (n,1) or a (1,n) vector
addParameter(ip, 'T', [], @(T) isnumeric(T) && (size(T,1)==1 || size(T,2)==1))

parse(ip, Q, t, varargin{:})
P = ip.Results.P;
PET = ip.Results.PET;
T = ip.Results.T;

% default setting reads as good data
error_flag = 0;
error_str = '';

% timestep checks
if isnumeric(t)
    error_flag = 1;
    t = datetime(t,'ConvertFrom','datenum');
    error_str = ['Warning: Converted datenum to datetime. ', error_str];
end

timesteps = diff(t);
timestep = median(timesteps);
if any(diff(timesteps)~=0) 
    error_flag = 1;
    error_str = ['Warning: Record is not continuous (some timesteps are missing). ', error_str];
end

% data checks
if min(Q)<0
    error_flag = 2;
    error_str = ['Error: Negative values in flow series. ', error_str];
    return
end

if all(Q==0)
    error_flag = 2;
    error_str = ['Error: Only zero flow in flow series. ', error_str];
    return
end

if length(Q) ~= length(t)
    error_flag = 2;
    error_str = ['Error: Flow series and time vector have different lengths. ', error_str];
    return
end

if any(isnan(Q))
    error_flag = 1;
    error_str = ['Warning: Ignoring NaNs in streamflow data. ', error_str];
end

if all(isnan(Q))
    error_flag = 2;
    error_str = ['Error: Only NaNs in streamflow data. ', error_str];
    return
end

if length(Q) < 30
    error_flag = 1;
    error_str = ['Warning: Extremely short time series. ', error_str];
end

% optionally check P
if ~isempty(P)
    
    if any(isnan(P))
        error_flag = 1;
        error_str = ['Warning: Ignoring NaNs in precipitation data. ', error_str];
    end
    
    if all(isnan(P))
        error_flag = 2;
        error_str = ['Error: Only NaNs in precipitation data. ', error_str];
        return
    end
    
    if length(Q) ~= length(P)
        error_flag = 2;
        error_str = ['Error: Precipitation and flow series have different lengths. ', error_str];
        return
    end
    
    if min(P)<0
        error_flag = 2;
        error_str = ['Error: Negative values in precipitation series. ', error_str];
        return
    end
    
end

% optionally check PET
if ~isempty(PET)
    
    if any(isnan(PET))
        error_flag = 1;
        error_str = ['Warning: Ignoring NaNs in potential evpotranspiration data. ', error_str];
    end
    
    if all(isnan(PET))
        error_flag = 2;
        error_str = ['Error: Only NaNs in potential evpotranspiration data. ', error_str];
        return
    end
    
    if length(Q) ~= length(PET)
        error_flag = 2;
        error_str = ['Error: Potential evpotranspiration and flow series have different lengths. ', error_str];
        return
    end
    
    if min(PET)<0
        error_flag = 1;
        error_str = ['Warning: Negative values in potential evpotranspiration series. ', error_str];
    end
    
end

% optionally check T
if ~isempty(T)
    
    if any(isnan(T))
        error_flag = 1;
        error_str = ['Warning: Ignoring NaNs in temperature data. ', error_str];
    end
    
    if all(isnan(T))
        error_flag = 2;
        error_str = ['Error: Only NaNs in temperature data. ', error_str];
        return
    end
    
    if length(Q) ~= length(T)
        error_flag = 2;
        error_str = ['Error: Temperature and flow series have different lengths. ', error_str];
        return
    end
    
    if min(T) < -273.15
        error_flag = 2;
        error_str = ['Error: Temperature cannot be less than the absolute minimum of -273.15 degC. ', error_str];
        return
    end
    
end

end
