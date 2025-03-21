function [A, phi, error_flag, error_str, fig_handles] = ...
    sig_SeasonalTranslation(Q, t, P, PET, varargin)
%sig_SeasonalTranslation calculates seasonal translation signatures.
%   Calculates amplitude ratio and phase shift between P-PET and Q. Uses
%   linear equation system to fit sine curve (see Gnann et al., 2020).
%
%   Notes:
%   The signatures were only properly tested in energy-limited catchments,
%   where AET is approximately equal to PET. AET can be estimated using a
%   simple water balance model, but that has not been properly tested yet.
%
%   INPUT
%   Q: streamflow [mm/timestep]
%   t: time [Matlab datetime]
%   P: precipitation [mm/timestep]
%   PET: potential evapotranspiration [mm/timestep]
%   OPTIONAL
%   period: period of cycle to be extracted [timestep], default = 365
%   useAET: whether to use AET in the input, default = false
%   plot_results: whether to plot results, default = false
%
%   OUTPUT
%   A: amplitude ratio [-]
%   phi: phase shift [timestep]
%   error_flag: 0 (no error), 1 (warning), 2 (error in data check), 3
%       (error in signature calculation)
%   error_str: string contraining error description
%   fig_handles: figure handles to manipulate figures (empty if plotting is
%       not requested)
%
%   EXAMPLE
%   % load example data
%   data = load('example/example_data/33029_daily.mat');
%   Q = data.Q;
%   t = data.t;
%   P = data.P;
%   PET = data.PET;
%   [A,phi] = sig_SeasonalTranslation(Q, t, P, PET, 'plot_results', true);
%
%   References
%   Gnann, S.J., Howden, N.J. and Woods, R.A., 2020. Hydrological
%   signatures describing the translation of climate seasonality into
%   streamflow seasonality. Hydrology and Earth System Sciences, 24(2),
%   pp.561-580.
%
%   Copyright (C) 2020
%   This software is distributed under the GNU Public License Version 3.
%   See <https://www.gnu.org/licenses/gpl-3.0.en.html> for details.

% check input parameters
if nargin < 4
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
% time series have to be numeric and either a (n,1) or a (1,n) vector
addRequired(ip, 'PET', @(PET) isnumeric(PET) && (size(PET,1)==1 || size(PET,2)==1))

% optional input arguments
addParameter(ip, 'period', 365, @isnumeric) % period of cycle (default 1 year)
addParameter(ip, 'useAET', false, @islogical) % whether to use AET in the input
addParameter(ip, 'plot_results', false, @islogical) % whether to plot results

parse(ip, Q, t, P, PET, varargin{:})
period = ip.Results.period;
useAET = ip.Results.useAET;
plot_results = ip.Results.plot_results;

% create empty figure handle
fig_handles = [];

% data checks
[error_flag, error_str, timestep, t] = util_DataCheck(Q, t, 'P', P, 'PET', PET);
if error_flag == 2
    A = NaN;
    phi = NaN;
    return
end

% get rid of NaN values (temporarily)
isn = (isnan(Q) | isnan(P) | isnan(PET)); % store NaN indices
% replace NaN days with mean to have a roughly closed water balance
Q(isn) = mean(Q,'omitnan');
P(isn) = mean(P,'omitnan');
PET(isn) = mean(PET,'omitnan');

% calculate signature
if useAET
    [~, AET] = util_StorageAndAET(Q, t, P, PET);
    F = P-AET; % proxy for incoming signal (forcing)
else
    F = P-PET; % proxy for incoming signal (forcing)
end
w = 2*pi/period; % angular frequency

% linear regression fitting
[A_F,phi_F,mean_F] = util_FitSineCurve(datenum(t), F, w);
[A_Q,phi_Q,mean_Q] = util_FitSineCurve(datenum(t), Q, w);

% calculate amplitude ratio
A = A_Q/A_F;

% calculate phase shift (we assume it's always less than 360 degrees)
phase_shift = phi_F-phi_Q;
if phase_shift < 0 % if Q "leads" P
    phase_shift = 2*pi+phase_shift;
end
phi = phase_shift./w;

if plot_results
    % get estimated sine curves
    F_hat = A_F*sin(w*datenum(t) + phi_F) + mean_F;
    Q_hat = A_Q*sin(w*datenum(t) + phi_Q) + mean_Q;
    fig = figure('Position',[100 100 700 250]); hold on
    % t_avg = datetime([0.*year(t),month(t),day(t)]);
    % t_avg(t_avg == datetime([0,12,31])) = NaT;
    plot(t,movmean(F,30),'color',[.8 .8 .8]);
    plot(t,movmean(Q,30),'color',[.4 .4 .4]);
    plot(t,F_hat,'-','Linewidth',1.5)
    plot(t,Q_hat,'-','Linewidth',1.5)
    xlabel('Date')
    ylabel('Flow [mm/timestep]')
    title(strcat('A [-]: ',num2str(round(A,2)),'; phi [timestep]: ',num2str(round(phi))))
    legend('Forcing','Streamflow','Sine curve forcing','Sine curve streamflow')
    fig_handles.SeasonalTranslation = fig;
end

end
