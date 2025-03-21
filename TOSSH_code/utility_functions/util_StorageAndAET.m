function [S, AET, fig_handles] = util_StorageAndAET(Q, t, P, PET, varargin)
%util_StorageAndAET calculates storage and actual evapotranspiration.
%   Calculates storage and actual evapotranspiration using a simple soil
%   moisture model (see e.g. Pfister et al., 2017).
%
%   Notes:
%   Includes a simple warm up to estimate field capacity and starts with
%   storage = field capacity. As an alternative, actual evapotranspiration
%   could be estimated separately as e.g. in Wlostowski et al. (2021).
%
%   INPUT
%   Q: streamflow [mm/timestep]
%   t: time [Matlab datenum]
%   P: precipitation [mm/timestep]
%   PET: potential evapotranspiration [mm/timestep]
%   OPTIONAL
%   field_capacity: field capacity [mm]
%
%   OUTPUT
%   S: storage [mm]
%   AET: actual evapotranspiration [mm/timestep]
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
%   [S, AET] = util_StorageAndAET(Q, t, P, PET);
%   [S, AET] = util_StorageAndAET(Q, t, P, PET, 'field_capacity', 100);
%
%   References
%   Pfister, L., Martínez-Carreras, N., Hissler, C., Klaus, J., Carrer,
%   G.E., Stewart, M.K. and McDonnell, J.J., 2017. Bedrock geology controls
%   on catchment storage, mixing, and release: A comparative analysis of 16
%   nested catchments. Hydrological Processes, 31(10), pp.1828-1845.
%   Wlostowski, A. N., Molotch, N., Anderson, S. P., Brantley, S. L.,
%   Chorover, J., Dralle, D., ... & Harman, C. (2021). Signatures of
%   hydrologic function across the critical zone observatory network.
%   Water Resources Research, 57(3), e2019WR026635.
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
addParameter(ip, 'field_capacity', [], @isnumeric) % field capacity for scaling PET to AET
addParameter(ip, 'plot_results', false, @islogical) % whether to plot results

parse(ip, Q, t, P, PET, varargin{:})
field_capacity = ip.Results.field_capacity;
plot_results = ip.Results.plot_results;

% create empty figure handle
fig_handles = [];

% get rid of NaN values (temporarily)
isn = (isnan(Q) | isnan(P) | isnan(PET)); % store NaN indices
Q(isn) = mean(Q,'omitnan');
P(isn) = mean(P,'omitnan');
PET(isn) = mean(PET,'omitnan');

% calculate storage and AET
AET = NaN(size(Q)); % actual evapotranspiration
S = NaN(size(Q)); % storage
a = NaN(size(Q)); % actual evapotranspiration scaling parameter

% if field capacity is not provided, we first need to estimate it
if isempty(field_capacity)

    field_capacity = 200; % field capacity initial guess
    AET(1) = PET(1);
    S(1) = field_capacity; % start at field capacity
    a(1) = 1;

    % warm up period to get estimate for initial S and field capacity
    for i = 2:length(t)
        if S(i-1) < field_capacity
            a(i) = S(i-1)/field_capacity;
        else
            a(i) = 1;
        end
        AET(i) = a(i)*PET(i);
        S(i) = P(i) - Q(i) - AET(i) + S(i-1);
        if S(i) < 0
            AET(i) = P(i) - Q(i) + S(i-1);
            S(i) = 0;
        end
        if AET(i) < 0
            AET(i) = 0;
        end
    end

    % use max as new field capacity
    field_capacity = max(S);
end

% begin with storage at field capacity
S(1) = field_capacity;

% now run loop again with provided or estimated field capacity
for i = 2:length(t)
    if S(i-1) < field_capacity
        a(i) = S(i-1)/field_capacity;
    else
        a(i) = 1;
    end
    AET(i) = a(i)*PET(i);
    S(i) = P(i) - Q(i) - AET(i) + S(i-1);
    if S(i) < 0
        AET(i) = P(i) - Q(i) + S(i-1);
        S(i) = 0;
    end
    if AET(i) < 0
        AET(i) = 0;
    end
end

% set Q etc back to NaN
Q(isn) = NaN;
P(isn) = NaN;
PET(isn) = NaN;
S(isn) = NaN;
AET(isn) = NaN;

% optional plotting
if plot_results
    fig = figure('Position',[100 100 700 250]); hold on;
    title('Storage and AET')
    h0 = plot(t, P, 'Color', [0.5 0.5 0.5]);
    h1 = plot(t, PET, 'c-');
    h2 = plot(t, AET, 'g-');
    h3 = plot(t, Q, 'b-');
    ylabel('Flux [mm/timestep]')
    yyaxis right
    ax = gca;
    ax.YAxis(2).Color = 'k';
    h4 = plot(t, S, 'r-');
    ylabel('Storage [mm]')
    legend([h0 h1 h2 h3 h4],{'Precipitation','PET','AET','Q','Storage'})
    xlabel('Date')
    fig_handles.StorageAndAET = fig;
end

end
