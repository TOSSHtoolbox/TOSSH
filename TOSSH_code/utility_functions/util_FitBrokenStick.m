function [err, fittedlines, slopes] = ...
    util_FitBrokenStick(breakpoints, x, y, zero_intercept)
%util_FitBrokenStick fits two or three part broken stick fit.
%   Fits two or three part broken stick fit (segmented regression) with
%   known break points. Returns norm of residuals to be used for
%   optimisation.
%   Adapted from online example by John D'Errico.
%
%   INPUT
%   breakpoints: breakpoints of broken stick (1 or 2 breakpoints)
%   x: independent variable
%   y: dependent variable
%   OPTIONAL
%   zero_intercept: should the intercept be 0?, default = false
%
%   OUTPUT
%   err: norm of residuals
%   fittedlines: start points and break points of fitted lines
%   slopes: slopes of the linear segments of the broken stick
%
%   EXAMPLE
%   x = [1:10]';
%   y = zeros(size(x));
%   y(1:5) = 2*x(1:5);
%   y(6:10) = y(5) + 0.5*(x(6:10)-5);
%   breakpoints = 5;
%   [err,fittedlines,slopes] = util_FitBrokenStick(breakpoints,x,y)
%
%   Copyright (C) 2020
%   This software is distributed under the GNU Public License Version 3.
%   See <https://www.gnu.org/licenses/gpl-3.0.en.html> for details.

if nargin < 3
    error('Not enough input arguments.')
end
if nargin < 4
    zero_intercept = false;
end

breakpoints = breakpoints(:).'; % make breakpoints into row vector

breaks = [min(x),sort(breakpoints),max(x)];
nx = length(x);
% which points lie in which interval?
xbins = discretize(x,breaks);
% breakpoints cannot be in the same interval
if numel(unique(xbins)) < numel(breakpoints) + 1
    err = 10^6;
    fittedlines = NaN(length(breakpoints)+2,2);
    slopes = NaN(length(breakpoints)+2,1);
    return
end

% write the problem in matrix form
if zero_intercept % intercept is zero
    if numel(breakpoints)==1
        A = [x - breaks(1),(x - breaks(2)).*(xbins == 2)];
    elseif   numel (breakpoints)==2
        A = [x - breaks(1),(x - breaks(2)).*(or(xbins == 2,xbins == 3)),(x - breaks(3)).*(xbins == 3)];
    else
        error('Function brokenstick only works for 1 or 2 breakpoints.')
    end
else % intercept will be optimised as well
    if numel(breakpoints)==1
        A = [ones(nx,1),x - breaks(1),(x - breaks(2)).*(xbins == 2)];
    elseif   numel (breakpoints)==2
        A = [ones(nx,1),x - breaks(1),(x - breaks(2)).*(or(xbins == 2,xbins == 3)),(x - breaks(3)).*(xbins == 3)];
    else
        error('Function brokenstick only works for 1 or 2 breakpoints.')
    end
end

% solve matrix
coef = A\y;
err = norm(y - A*coef);
if zero_intercept
    coef = [0; coef];
end

% unpack the coefficients
c1 = coef(1);
s = [coef(2:end)];
css = cumsum(s);
s = [s;s(end)];
css = [css;css(end)];
bi = [breaks(1:end-1).';breaks(end-1)];
fittedlines = zeros(length(breaks),2);
slopes = zeros(length(breaks),1);
fittedlines(:,1) = breaks.';
fittedlines(1:2,2) = (c1 - breaks(1).*s(1)) + s(1).*breaks(1:2).';
fittedlines(3,2) = (c1 - breaks(1).*s(1) - breaks(2).*s(2)) + (s(1)+s(2)).*breaks(3).';
if numel (breakpoints)==2
    fittedlines(4,2) = (c1 - sum(breaks(1:3).'.*s(1:3))) + sum(s(1:3)).*breaks(4).';
end
slopes = cumsum(s);slopes=slopes(1:end-1);

end


