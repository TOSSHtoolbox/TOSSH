function [RecessionParameters_a, RecessionParameters_b, RecessionParameters_T0] = ...
    util_Recession_PostProcess(RecessionParametersTemp, Q_cell, i, varargin)
%util_Recession_PostProcess post-processes recession parameters.
%
%   Computes representative recession parameters (a, b, T0) from temporary
%   estimates using median-based approaches.
%
%   RecessionParameters_T0 is characteristic timescale of recessions at median flow; It can be obtained
%   by fitting a line to the dQ/dt versus Q point cloud in log-log space for each individual recession,
%   with Q scaled by median Q; T0 is the median value of −1/intercept (McMillan et al., 2021).
%
%   RecessionParameters_T0 can be derived from RecessionParameters_a and RecessionParameters_b as follows (McMillan et al, 2014):
%   Let Qhat = Q/Qmedian, then:
%       dQhat/dt = - (1/T0) * Qhat^b
%   Separating constant terms,
%       1/Qmedian * dQ/dt = - (1/T0) * Q^b * (1/Qmedian)^b
%       dQ/dt = - (1/T0) * Q^b * (1/Qmedian) ^ (b-1)
%   As our estimate of a and b parameters are from:
%       dQ/dt = - a * Q^b
%   It leads to:
%       a = (1/T0) * (1/Qmedian) ^ (b-1)
%       T0 = (1/a) * (1/Qmedian) ^ (b-1)
%       T0 = 1 / (a * Qmedian ^ (b-1)) 
%
%   INPUT
%   RecessionParametersTemp: [n x 2] matrix with columns [a, b]
%   Q_cell: cell array with discharge values
%   i: index for current time series
%
%   OPTIONAL
%   'useMedianPair': logical, default = false
%       true  -> use median b and corresponding (a,b) pair
%       false -> use separate medians for a and b
%
%   OUTPUT
%   RecessionParameters_a  : representative recession coefficient a [-]
%   RecessionParameters_b  : representative recession exponent b [-]
%   RecessionParameters_T0: characteristic recession timescale at median discharge [timestep]
%
%   References
%   McMillan, H., Gueguen, M., Grimon, E., Woods, R., Clark, M., & Rupp, D. E. (2014).
%   Spatial variability of hydrological processes and model structure diagnostics in a 50 km2 catchment.
%   Hydrological Processes, 28(18), 4896-4913. https://doi.org/10.1002/hyp.9988
%
%   Copyright (C) 2026
%   This software is distributed under the GNU Public License Version 3.
%   See <https://www.gnu.org/licenses/gpl-3.0.en.html> for details.

if nargin < 3
    error('Not enough input arguments.')
end

ip = inputParser;
ip.CaseSensitive = true;

addRequired(ip, 'RecessionParametersTemp', @isnumeric)
addRequired(ip, 'Q_cell')
addRequired(ip, 'i', @isnumeric)

addParameter(ip, 'useMedianPair', false, @islogical)

parse(ip, RecessionParametersTemp, Q_cell, i, varargin{:})
useMedianPair = ip.Results.useMedianPair;

% median discharge 
Qmed = median(Q_cell{i}(Q_cell{i} > 0), 'omitnan');

% parameter selection 
if useMedianPair
    % median b and corresponding (a,b) pair
    b_med = median(RecessionParametersTemp(:,2), 'omitnan');
    [~, ind] = min(abs(RecessionParametersTemp(:,2) - b_med));

    RecessionParameters_a = RecessionParametersTemp(ind,1);
    RecessionParameters_b = RecessionParametersTemp(ind,2);

    % T0 computed from selected pair
    RecessionParameters_T0 = 1 ./ ...
        (RecessionParameters_a .* Qmed.^(RecessionParameters_b - 1));

else
    % separate medians
    RecessionParameters_a = median(RecessionParametersTemp(:,1), 'omitnan');
    RecessionParameters_b = median(RecessionParametersTemp(:,2), 'omitnan');

    % T0 from distribution
    T0_temp = 1 ./ ...
        (RecessionParametersTemp(:,1) .* Qmed.^(RecessionParametersTemp(:,2) - 1));

    % filter for reasonable b values
    idxValid = (RecessionParametersTemp(:,2) > 0.5) & ...
               (RecessionParametersTemp(:,2) < 5);

    RecessionParameters_T0 = median(T0_temp(idxValid), 'omitnan');
end

end