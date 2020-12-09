function [X_sub, t_sub] = util_ExtractSubPeriod(X, t, subperiod, option)
%util_ExtractSubPeriod extracts subperiod (months) from time series.
%   Extracts certain months from time series, e.g. the low flow season and 
%   the high flow season. See also Euser et al. (2013) who use May to
%   September as low flow season and November to April as high flow season.
%
%   INPUT
%   X: time series - should be continuous (e.g. Q, typically [mm/timestep])
%   t: time [Matlab datetime]
%   subperiod: subperiod to be extracted, e.g. [11:12, 1:4] for November
%       to April or [5:9] for May to September
%   OPTIONAL
%   option: option to delete values outside subperiod ('delete') or assign
%       NaN to values outside subperiod ('nan')
%
%   OUTPUT
%   X_sub: subperiod time series
%   t_sub: subperiod time
%
%   EXAMPLE
%   % load example data
%   data = load('example/example_data/33029_daily.mat');
%   Q = data.Q;
%   t = data.t;
%   subperiod = [11:12, 1:4];
%   [Q_sub, t_sub] = ...
%       util_ExtractSubPeriod(Q, t, subperiod);
%
%   References
%   Euser, T., Winsemius, H.C., Hrachowitz, M., Fenicia, F., Uhlenbrook, S.
%   and Savenije, H.H.G., 2013. A framework to assess the realism of model
%   structures using hydrological signatures. Hydrology and Earth System 
%   Sciences, 17 (5), 2013.
%   
%   Copyright (C) 2020
%   This software is distributed under the GNU Public License Version 3.
%   See <https://www.gnu.org/licenses/gpl-3.0.en.html> for details.

if nargin < 4
    option = 'delete';
end

% get months
[~, month_vec, ~] = ymd(t);

if strcmp(option,'delete')
    % option 1 - delete values outside subperiod
    t_sub = t(ismember(month_vec,subperiod));
    X_sub = X(ismember(month_vec,subperiod));
elseif strcmp(option,'nan')
    % option 2 - assign NaN to values outside subperiod
    t_sub = t;
    X_sub = X;
    X_sub(ismember(month_vec,subperiod)) = NaN;
else
    error('Subperiod option specified incorrectly, the options are delete or nan.')
end

%{
figure; hold on;
plot(t,X);
plot(t_sub,X_sub,'--')
xlabel('Date')
ylabel('Flow')
%}

end