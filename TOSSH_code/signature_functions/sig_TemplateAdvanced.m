function [ExampleSignature, error_flag, error_str, fig_handles] = ...
    sig_TemplateAdvanced(Q, t, P, PET, T, param, varargin)
%sig_TemplateAdvanced calculates [Enter brief description of the signature].
%   [Enter a more detailed description of the signature, possibly including
%   relevant references and information about different options.]
%
%   INPUT
%   Q: streamflow [mm/timestep]
%   t: time [Matlab datetime]
%   P: precipitation [mm/timestep]
%   PET: potential evapotranspiration [mm/timestep]
%   T: temperature [�C]
%   param: required parameter
%   [...]
%   OPTIONAL
%   opt_param: optional parameter
%   plot_results: whether to plot results, default = 0
%   [...]
%
%   OUTPUT
%   example_signature: example signature [-]
%   error_flag: 0 (no error), 1 (warning), 2 (error in data check), 3
%       (error in signature calculation)
%   error_str: string contraining error description
%   fig_handles: figure handles to manipulate figures (empty if plotting is
%       not requested)
%   [...]
%   
%   EXAMPLE
%   % load example data
%   data = load('example/example_data/33029_daily.mat');
%   Q = data.Q;
%   t = data.t;
%   P = data.P;
%   PET = data.PET;
%   T = data.T;
%   param = 1;
%   ExampleSignature = sig_TemplateAdvanced(Q,t,P,PET,T,param);
%   [...]
%
%   References
%   [...]
%
%   Copyright (C) 2020
%   This software is distributed under the GNU Public License Version 3.
%   See <https://www.gnu.org/licenses/gpl-3.0.en.html> for details.

% check input parameters
% [Change if there is a different number of required inputs.]
if nargin < 6
    error('Not enough input arguments.')
end

ip = inputParser;
ip.CaseSensitive = true; 

% required input arguments
%[Delete the input parameters you don't need.]
% time series have to be numeric and either a (n,1) or a (1,n) vector
addRequired(ip, 'Q', @(Q) isnumeric(Q) && (size(Q,1)==1 || size(Q,2)==1))
% date time series has to be numeric or datetime and either a (n,1) or a (1,n) vector
addRequired(ip, 't', @(t) (isnumeric(t) || isdatetime(t)) && (size(t,1)==1 || size(t,2)==1))
% time series have to be numeric and either a (n,1) or a (1,n) vector
addRequired(ip, 'P', @(P) isnumeric(P) && (size(P,1)==1 || size(P,2)==1))
% time series have to be numeric and either a (n,1) or a (1,n) vector
addRequired(ip, 'PET', @(PET) isnumeric(PET) && (size(PET,1)==1 || size(PET,2)==1))
% time series have to be numeric and either a (n,1) or a (1,n) vector
addRequired(ip, 'T', @(T) isnumeric(T) && (size(T,1)==1 || size(T,2)==1))
% param has to be numeric and either a (n,1) or a (1,n) vector
addRequired(ip, 'param', @(param) isnumeric(param) && (size(param,1)==1 || size(param,2)==1))

% optional input arguments
addParameter(ip, 'opt_param', false, @islogical) % optional parameter
addParameter(ip, 'plot_results', false, @islogical) % whether to plot results

parse(ip, Q, t, P, PET, T, param, varargin{:}) %[Delete the input parameters you don't need.]
opt_param = ip.Results.opt_param;
plot_results = ip.Results.plot_results;

% create empty figure handle
fig_handles = [];

% data checks
[error_flag, error_str, timestep, t] = util_DataCheck(Q, t, 'P', P, 'PET', PET, 'T', T); %[Delete the input parameters you don't need.]
if error_flag == 2
    ExampleSignature = NaN;
    return
end

% [Add warnings/errors to indicate potentially problematic inputs here.]
if opt_param
    error_flag = 1;
    error_str = ['Warning: You have set the optional parameter to true. ', error_str];
end

% calculate signature
ExampleSignature = param;
% [Add well commented signature here.]
% ...
% ...

% optional plotting
if plot_results
    fig = figure('pos',[100 100 350 300]); hold on
    plot(ExampleSignature);
    fig_handles.TemplateAdvanced = fig;
end

end
