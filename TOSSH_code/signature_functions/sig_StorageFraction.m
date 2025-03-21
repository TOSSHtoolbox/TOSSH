function [S_fraction, S_active, S_total, error_flag, error_str, fig_handles] = ...
    sig_StorageFraction(Q, t, P, PET, varargin)
%sig_StorageFraction calculates ratio between active and total storage.
%   Maximum storage deficit in the series (active storage) and total
%   storage volume (extrapolation to find storage deficit at near-zero
%   flow; see Pfister et al., 2017).
%
%   Notes:
%   This is our own implementation of Pfister et al. (2017) and it has not
%   been tested for a wide range of catchments. Please investigate the plot
%   to see if the signatures are calculated reasonably.
%
%   INPUT
%   Q: streamflow [mm/timestep]
%   t: time [Matlab datetime]
%   P: precipitation [mm/timestep]
%   PET: potential evapotranspiration [mm/timestep]
%   OPTIONAL
%   field_capacity: field capacity [mm]
%   bin_size: bin size used to determine envelope, default: 10 mm
%   fit_range: range of envelope to which linear line should be fitted,
%              default: 50th to 99th percentile, i.e. [0.5 0.99]
%   plot_results: whether to plot results, default = false
%
%   OUTPUT
%   S_fraction: ratio between active and total storage capacity [-]
%   S_active: active storage capacity [mm]
%   S_total: total storage capacity [mm]
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
%   S_fraction = sig_StorageFraction(Q, t, P, PET);
%   [S_fraction, S_active, S_total] = sig_StorageFraction(Q, t, P, PET, 'plot_results', true);
%
%   References
%   Pfister, L., Martínez-Carreras, N., Hissler, C., Klaus, J., Carrer,
%   G.E., Stewart, M.K. and McDonnell, J.J., 2017. Bedrock geology controls
%   on catchment storage, mixing, and release: A comparative analysis of 16
%   nested catchments. Hydrological Processes, 31(10), pp.1828-1845.
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
addParameter(ip, 'field_capacity', [], @isnumeric) % field capacity
addParameter(ip, 'bin_size', 10, @isnumeric) % bin size used to determine envelope
addParameter(ip, 'fit_range', [0.5 0.99], @isnumeric) % range of envelope to which linear line should be fitted
addParameter(ip, 'plot_results', false, @islogical) % whether to plot results

parse(ip, Q, t, P, PET, varargin{:})
field_capacity = ip.Results.field_capacity;
bin_size = ip.Results.bin_size;
fit_range = ip.Results.fit_range;
plot_results = ip.Results.plot_results;

% create empty figure handle
fig_handles = [];

% data checks
[error_flag, error_str, timestep, t] = util_DataCheck(Q, t, 'P', P, 'PET', PET);
if error_flag == 2
    S_fraction = NaN;
    S_active = NaN;
    S_total = NaN;
    return
end

% calculate signature

% !NaN check now in util function!
% % get rid of NaN values (temporarily)
% isn = (isnan(Q) | isnan(P) | isnan(PET)); % store NaN indices
% % replace NaN days with mean to have a roughly closed water balance
% Q(isn) = mean(Q,'omitnan');
% P(isn) = mean(P,'omitnan');
% PET(isn) = mean(PET,'omitnan');

% estimate storage
[S, ~] = util_StorageAndAET(Q, t, P, PET, 'field_capacity', field_capacity);

% Q(isn) = NaN; % set Q and corresponding S back to NaN
% S(isn) = NaN;

S_max = max(S); % maximum storage capacity
D = S_max - S; % storage deficit

% get FDC and corresponding deficits
[Q_sorted,ind] = sort(Q);
% Q_ranked = tiedrank(Q_sorted);
Q_ranked = [1:length(Q)]'; % give unique (random) rank to every measurement
FDC = 1 - Q_ranked./length(Q_ranked); % flow duration curve with unique ranks
D_sorted = D(ind);

% get deficit value determined from the 99th percentile of the observed FDC
indices = 1:length(FDC);
ind99 = max(indices(FDC >= 0.99));
S_active = D_sorted(ind99);

% calculate total storage capacity S_total
D_rounded = round((1/bin_size).*D_sorted).*bin_size; % round deficit and get unique values
D_rounded(isnan(D_rounded)) = []; % remove NaN values
D_unique = unique(D_rounded);
D_envelope = NaN(size(D_unique));
Q_envelope = NaN(size(D_unique));
for i = 1:length(D_unique)
    Q_tmp = Q_sorted(D_rounded==D_unique(i));
    D_tmp = D_sorted(D_rounded==D_unique(i));
    [Q_envelope(i),index] = min(Q_tmp); % get min of bin to obtain envelope points
    D_envelope(i) = D_tmp(index); % get corresponding deficits
end

% fit line to envelope (default: use 50th to 100th percentile (upper limit))
low = round(length(Q_envelope)*fit_range(1));
upp = round(length(Q_envelope)*fit_range(2));

% ephemeral catchments need special treatment
if min(Q) < 0.001
    S_total = max(D);
    S_active = S_total;
    error_flag = 1;
    error_str = ['Warning: Storage ratio calculation unreliable for ephemeral catchments. ', error_str];
else
    [m, c, ~] = util_FitLinear(Q_envelope(low:upp),D_envelope(low:upp));
    % todo: add uncertainty
    S_total = m + c.*0.001;
    if c >= 0 % slope of fitted line has to be negative
        c = NaN;
        S_total = NaN;
        error_flag = 1;
        error_str = ['Warning: Total storage could not be estimated properly. ', error_str];
    end
end

% check results
if isempty(S_active) || isempty(S_total)
    S_fraction = NaN;
    S_total = NaN;
    S_active = NaN;
    error_flag = 3;
    error_str = ['Error: Active or total storage could not be calculated. ', error_str];
    return
end

if S_active > S_total
    error_flag = 1;
    error_str = ['Warning: Estimated active storage is larger than total storage. ', error_str];
end

S_fraction = S_active/S_total; % calculate storage fraction

% optional plotting
if plot_results
    fig = figure('Position',[100 100 350 300]); hold on
    p1 = plot(Q,D,'-','color',[.9 .9 .9]);
    p2 = plot(Q,D,'.','color',[.5 .5 .5]);
    p3 = plot(Q_envelope,D_envelope,'k o','linewidth',2);
    if min(Q) < 0.001
        p4 = plot(NaN,NaN);
    else
        p4 = plot([median(Q_envelope):-0.001:0.001],...
            m+[median(Q_envelope):-0.001:0.001].*c,'g-','linewidth',1.5);
    end
    p5 = plot(0.001,S_total,'g x','linewidth',1.5);
    p6 = plot(Q_sorted(ind99),S_active,'r o','linewidth',2);
    xlabel('Flow [mm/d]')
    ylabel('Storage deficit [mm]')
    legend([p2 p3 p5 p6 p4],...
        {'All points','Envelope','S_{total}','S_{active}','Fitted Line'})
    %set(gca,'xScale','log');
    fig_handles.StorageFraction = fig;
end

end
