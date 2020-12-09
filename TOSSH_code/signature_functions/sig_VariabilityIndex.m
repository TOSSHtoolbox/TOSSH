function [VariabilityIndex, error_flag, error_str] = sig_VariabilityIndex(Q, t)
%sig_VariabilityIndex calculates variability index (VI) from FDC.
%   VI is the standard deviation of the common logarithms of discharge 
%   determined at 10% intervals from 10% to 90% of the cumulative frequency 
%   distribution (flow duration curve, FDC). Low variability index shows 
%   higher water storage (Estrany et al., 2010).
%
%   INPUT
%   Q: streamflow [mm/timestep]
%   t: time [Matlab datetime]
%
%   OUTPUT
%   VariabilityIndex: variability index [-]
%   error_flag: 0 (no error), 1 (warning), 2 (error in data check), 3
%       (error in signature calculation)
%   error_str: string contraining error description
%
%   EXAMPLE
%   % load example data
%   data = load('example/example_data/33029_daily.mat');
%   Q = data.Q;
%   t = data.t;
%   VariabilityIndex = sig_VariabilityIndex(Q,t);
%
%   References
%   Estrany, J., Garcia, C. and Batalla, R.J., 2010. Hydrological response
%   of a small mediterranean agricultural catchment. Journal of Hydrology,
%   380(1-2), pp.180-190.
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

parse(ip, Q, t)

% data checks
[error_flag, error_str, timestep, t] = util_DataCheck(Q, t);
if error_flag == 2
    VariabilityIndex = NaN;
    return
end

% calculate signature
% get ranks as a proxy for exceedance probabilities
Q_sorted = sort(Q,'descend');

% percentiles required are 10%, 20%, ..., 90%
percs = [10:10:90];

% get the corresponding rank of the FDC values
indices_percs = round(length(Q_sorted).*percs./100);

% get the flow value at each rank
flow_percs = Q_sorted(indices_percs);

% variation needed if some flow percentiles are zero for an intermittent
% stream - exclude these from the calculation
recs = flow_percs > 0;

% VI is the standard deviation of the common logarithms of discharge 
% determined at 10% intervals from 10% to 90% of the cumulative frequency 
% distribution
VariabilityIndex = std(log10(flow_percs(recs)));

end

