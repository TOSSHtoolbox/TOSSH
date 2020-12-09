function [dQdt, Qm, flow_section, R2] = util_dQdt(Q, t, flow_section, varargin)
%util_dQdt calculates flow rate gradient with various options.
%   Default method is the "classical" Brutsaert and Nieber (1977) method.
%   The Roques et al. (2017) method (ETS) is suggested for a more robust
%   calculation of dQ/dt, but it is not the default as it is rather slow.
%
%   INPUT
%   Q: streamflow [mm/timestep]
%   t: time [Matlab datenum]
%   flow_section: n-by-2 array where n is the number of recession segments
%   OPTIONAL
%   method: method for dQdt calculation
%
%   OUTPUT
%   dQdt: flow rate gradient [mm/timestep^2]
%   Qm: corresponding flow [mm/timestep]
%   flow_section: updated flow_section array (some recession points have to
%       be removed due to approx. of derivative)
%   R2: R^2 from exponential time stepping method
%
%   EXAMPLE
%   % load example data 
%   data = load('example/example_data/33029_daily.mat'); 
%   Q = data.Q; 
%   t = data.t;
%   flow_section = util_RecessionSegments(Q,t); % get recession segments
%   [dQdt, Qm, flow_section, R2] = util_dQdt(Q, t, flow_section);
%   [dQdt, Qm, flow_section, R2] = util_dQdt(Q, t, flow_section, 'method', 'ETS');
%
%   References
%	Brutsaert, W. and Nieber, J.L., 1977. Regionalized drought flow 
%   hydrographs from a mature glaciated plateau. Water Resources Research, 
%   13(3), pp.637-643.
%   Thomas, B.F., Vogel, R.M. and Famiglietti, J.S., 2015. Objective 
%   hydrograph baseflow recession analysis. Journal of hydrology, 525, 
%   pp.102-112.
%   Roques, C., Rupp, D.E. and Selker, J.S., 2017. Improved streamflow 
%   recession parameter estimation with attention to calculation of dQ/dt.
%   Advances in Water Resources, 108, pp.29-43.
%
%   Copyright (C) 2020
%   This software is distributed under the GNU Public License Version 3.
%   See <https://www.gnu.org/licenses/gpl-3.0.en.html> for details.

% check input parameters
if nargin < 3
    error('Not enough input arguments.')
end

ip = inputParser;
ip.CaseSensitive = true; 

% required input arguments
% time series have to be numeric and either a (n,1) or a (1,n) vector
addRequired(ip, 'Q', @(Q) isnumeric(Q) && (size(Q,1)==1 || size(Q,2)==1)) 
% date time series has to be numeric or datetime and either a (n,1) or a (1,n) vector
addRequired(ip, 't', @(t) (isnumeric(t) || isdatetime(t)) && (size(t,1)==1 || size(t,2)==1)) 
addRequired(ip, 'flow_section', @(flow_section) isnumeric(flow_section) && size(flow_section,2)==2) % recession segments 

% optional input arguments
addParameter(ip, 'method', 'BN', @ischar) 

parse(ip, Q, t, flow_section, varargin{:})
method = ip.Results.method;

dQdt = NaN(size(Q));
Qm = NaN(size(Q));
m = NaN(size(Q));
R2 = ones(size(Q)); % weights

switch method    
        
    case 'BN' % Brutsaert and Nieber (1979)
        for j = 1:size(flow_section,1)
            rec = [flow_section(j,1):flow_section(j,2)]'; % get recession
            dQdt(rec(2:end)) = (Q(rec(2:end))-Q(rec(1:end-1)))./days(t(2)-t(1));
            Qm(rec(2:end)) = (Q(rec(2:end)) + Q(rec(1:end-1)))./2;            
            flow_section(j,1) = flow_section(j,1)+1; % shorten recession
        end
       
    case 'backwards' % similar to Brutsaert and Nieber (1979), but we keep measured Q, see also Thomas et al. (2015)
        for j = 1:size(flow_section,1)
            rec = [flow_section(j,1):flow_section(j,2)]'; % get recession
            dQdt(rec(2:end)) = (Q(rec(2:end))-Q(rec(1:end-1)))./days(t(2)-t(1));
            Qm(rec(2:end)) = Q(rec(2:end));
            flow_section(j,1) = flow_section(j,1)+1; % shorten recession
        end
        
    case 'ETS' % exponential time stepping following Roques et al. (2017)        
        for j = 1:size(flow_section,1)
            rec = [flow_section(j,1):flow_section(j,2)]'; % get recession
            n = 0.1*(length(rec)); % n = 10% of recession led to good results according to Roques et al. (2017) 
            gamma = util_FitExponential(Q(rec), t(rec)); % get gamma
            m(rec) = 1 + ceil(n.*exp(-1./(gamma.*[1:length(rec)]))); 
            
            i = rec(1);
            while i+m(i) <= rec(end)
                Qm(i) = mean(Q(i:i+m(i)));
                [~, dQdt(i), R2(i)] = util_FitLinear(...
                    datenum(t(i:i+m(i))),Q(i:i+m(i)));
                i = i+1;
            end 
            
            flow_section(j,2) = flow_section(j,2) - m(i);%(m(i)+1); % shorten recession            
        end
        
    otherwise
        error('Differentiation method not available.')
end

% dQdt has to be negative and weights have to be non-negative
dQdt(dQdt>=0 | isinf(dQdt)) = NaN;
Qm(isnan(dQdt)) = NaN;
R2(R2<=0) = 1e-18; % practically 0 weight
R2(isnan(R2)) = 1e-18;

end

