function [A, phi, k] = util_FitSineCurve(x, y, w)
%util_FitSineCurve fits sine curve to time series.
%   y = A*sin(w*x + phi) + k;
%
%   INPUT
%   x: independent variable (typically time)
%   y: dependent variable (typically flow)
%   w: angular frequency
%
%   OUTPUT
%   A: amplitude
%   phi: phase
%   k: offset
%
%   EXAMPLE
%   x = [0:5*365]';
%   w = 2*pi/365;
%   y = 1 + sin(w.*x + pi/2);
%   [A,phi,k] = util_FitSineCurve(x,y,w)
%
%   Copyright (C) 2020
%   This software is distributed under the GNU Public License Version 3.
%   See <https://www.gnu.org/licenses/gpl-3.0.en.html> for details.

% create matrix
M = ones(length(x),3);
M(:,2) = cos(w*x);
M(:,3) = sin(w*x);

% change NaN values to median for linear regression
% y(isnan(y)) = median(y,'omitnan');

% solve equation system
b = M\y;

% get estimated sine curve parameters
phi = atan2(b(2),b(3)); % get unambigous value for phi
A = sqrt(b(2)^2+b(3)^2);
k = b(1);

%{
% get estimated sine curve
y_hat = A*sin(w*x + phi) + k;
figure; 
plot(x,y); hold on
plot(x,y_hat,'r --')
%}

end