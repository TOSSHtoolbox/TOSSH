function [Q_b] = util_UKIH_Method(Q, varargin)
%util_UKIH_Method estimates baseflow with UKIH method.
%   Estimates baseflow UKIH "smoothed minima" method (UK Institute of
%   Hydrology, 1980).
%
%   INPUT
%   Q: [mm/timestep]
%   OPTIONAL
%   n_days: length of data blocks, default = 5 days
%   
%   OUTPUT
%   Q_b: baseflow [mm/timestep]
%
%   EXAMPLE
%   % load example data 
%   data = load('example/example_data/33029_daily.mat'); 
%   Q = data.Q; 
%   t = data.t; 
%   Q_b = util_UKIH_Method(Q);
%   Q_b90 = util_UKIH_Method(Q, 'n_days',  90);
%
%   References
%   UK Institute of Hydrology (Great Britain), 1980. Low Flow Studies
%   Reports. Institute of Hydrology.
%
%   Copyright (C) 2020
%   This software is distributed under the GNU Public License Version 3.
%   See <https://www.gnu.org/licenses/gpl-3.0.en.html> for details.

if nargin < 1
    error('Not enough input arguments.')
end

ip = inputParser;
ip.CaseSensitive = true; 

% required input arguments
% time series have to be numeric and either a (n,1) or a (1,n) vector
addRequired(ip, 'Q', @(Q) isnumeric(Q) && (size(Q,1)==1 || size(Q,2)==1))

% optional input arguments
addParameter(ip, 'n_days', 5, @isnumeric)

parse(ip, Q, varargin{:})
n_days = ip.Results.n_days;

if floor(n_days)~=n_days && n_days<1
    error('Filter window must be an integer larger than zero.')
end

% Baseflow separation is problematic with NaN values. Therefore, we set NaN
% values to median, apply the filter, and then set baseflow to NaN where
% streamflow is NaN. If there are a lot of NaN values, we encourage the
% user to either interpolate these values or to calculate the signature for
% each block individually and then calculate a weighted average.
Q_tmp = Q;
Q_tmp(isnan(Q)) = median(Q,'omitnan');

% calculate baseflow
[Q_b, t_ind] = UKIH_Method(Q_tmp, n_days); % 5 is the default parameter

% use minimum baseflow to fill in missing values at the beginning
B_tmp = min(Q_tmp)*ones(size(Q_tmp));
B_tmp(t_ind) = Q_b;
Q_b = B_tmp;

% set baseflow to NaN where streamflow is NaN
Q_b(isnan(Q)) = NaN;

end

function [Q_b, t_ind] = UKIH_Method(Q, n_days)
%UKIH_Method Helper function that runs the UKIH method.

n = length(Q);
Q_min5 = NaN(round(n/n_days),1); % 5 day minima
min_i = NaN(round(n/n_days),1); % corresponding indices
ind5 = 1; % minima counter
TP = 0; % turning points
t_TP = 0; % corresponding time/index
indTP = 1; % TP counter

for i = 1:n_days:floor(n/n_days)*n_days % divide in non-overlapping n-day blocks
    [Q_min5(ind5), min_i(ind5)] = min(Q(i:i+(n_days-1))); % find minimum
    if ind5 <= 2 % need at least three minima
    elseif      Q_min5(ind5-1)*0.9 < Q_min5(ind5-2) ...
            &&  Q_min5(ind5-1)*0.9 < Q_min5(ind5) % check if baseflow ordinate
        TP(indTP) = Q_min5(ind5-1);
        t_TP(indTP) = i - n_days - 1 + min_i(ind5-1); % get corresponding index
        indTP = indTP + 1;
    end
    ind5 = ind5 + 1;
end

t_ind = [t_TP(1):t_TP(end)]';
if t_ind == 0
    t_ind = [1:length(Q)]';
    Q_b = NaN(size(Q));
else
    Q_b = interp1q(t_TP',TP',t_ind); % linear interpolation
    Qt = Q(t_ind);
    Q_b(Q_b>Qt) = Qt(Q_b>Qt); % constrain B, so that B is never larger than Q
end

end

