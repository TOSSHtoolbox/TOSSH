function [a, b, R2] = util_FitLinear(x,y)
%util_FitLinear fits linear function and returns parameters and residuals.
%   y = a + b*x;
%   
%   INPUT
%   x: independent variable
%   y: dependent variable
%
%   OUTPUT
%   a: offset parameter
%   b: slope parameter
%   R2: residuals
%   
%   EXAMPLE
%   x = [0:5];
%   y = 1 + 2*x;
%   [a, b, R2] = util_FitLinear(x,y)
%
%   Copyright (C) 2020
%   This software is distributed under the GNU Public License Version 3.
%   See <https://www.gnu.org/licenses/gpl-3.0.en.html> for details.

n = length(x);
SSxy = sum(x.*y) - sum(x)*sum(y)/n;
SSxx = sum(x.^2) - sum(x)^2/n;
b = SSxy/SSxx;
a = mean(y) - b*mean(x);
y_hat = a + b*x;
R2 = 1 - sum((y - y_hat).^2)./sum((y - mean(y)).^2);

end
