function [thresh, slope, slope_linear, p_value] = util_Threshold(x, y, varargin)
%util_Threshold fits threshold to data and returns significance value.
%   This function fits a threshold to point data and returns a
%   significance value that indicates whether threshold exists. Different 
%   types of threshold shape are possible and should be specified in the 
%   shape parameter.
%
%   INPUT
%   x: x variable
%   y: y variable
%   OPTIONAL
%   shape: shape of threshold
%
%   OUTPUT
%   thresh: threshold
%   slope: slope of line starting at threshold
%   slope_linear: slope of linear function strating at origin
%   p_value: p_value < 0.05 implies that the threshold exists
%
%   EXAMPLE
%   x = [1:10]';
%   y = zeros(size(x));
%   y(1:5) = 0*x(1:5);
%   y(6:10) = y(5) + 0.5*(x(6:10)-5);
%   [thresh,slope,slope_linear,p_value] = util_Threshold(x, y);
%
%   Copyright (C) 2020
%   This software is distributed under the GNU Public License Version 3.
%   See <https://www.gnu.org/licenses/gpl-3.0.en.html> for details.

% check input parameters
if nargin < 2
    error('Not enough input arguments.')
end

ip = inputParser;
ip.CaseSensitive = true;

% required input arguments
% values to fit threshold for have to be numeric and either a (n,1) or a (1,n) vector
addRequired(ip, 'x', @(x) isnumeric(x) && (size(x,1)==1 || size(x,2)==1))
% date time series has to be numeric or datetime and either a (n,1) or a (1,n) vector
addRequired(ip, 'y', @(y) isnumeric(y) && (size(y,1)==1 || size(y,2)==1))

% optional input arguments
% shape of threshold. zero_linear fits a flat zero portion followed by a linear segment
addParameter(ip, 'shape', 'zero_linear', @ischar) 

parse(ip, x, y, varargin{:})
shape = ip.Results.shape;

% data checks
if length(x) ~= length(y)
    error('x and y data for threshold fit must have equal lengths')
end

% 1. compute goodness of fit for trial run of possible breakpoints
n = 30;
% calculate potential breakpoint locations equally spaced between maximum
% and minimum of data
bil = linspace(.5/n,1-.5/n,n);
bi = min(x)+bil.*(max(x)-min(x));
% calculate potential slopes equally spaced between 0 and max
% where max slope is for a line starting 80% of the way along the x axis 
% and reaching the highest point
max_slope = max(y)/(0.8*max(x));
si = linspace(max_slope/n,1-max_slope/n,n).*max_slope;

% possible breakpoint locations depend on shape of threshold
% create an array of possible breakpoint locations
switch shape
    case 'zero_linear'
        % test the broken stick fit for each possible breakpoint location
        trialfit = zeros(n,1);
        for i = 1:n
            trialfit(i) = FitBrokenLine([bi(i),si(i)],x,y);
        end
        [~,i] = min(trialfit(:)); % find the best of the breakpoint/slopes tested
        b0 = [bi(i) si(i)]; % initialize breakpoint/slope to that value
end

% create a helper function to pass the x, y data to the objective function
fun = @(b)FitBrokenLine(b,x,y);

% 2. search for the optimal breakpoint, starting from the best so far, b0
% b_old=fminsearch(fun,b0) % returns b, the optimal breakpoint/slope
options = optimoptions('fmincon','Display','off');
b = fmincon(fun,b0,[],[],[],[],[0, 0],[inf, inf],[],options);

% 3. call brokenline to find the sse for threshold vs non-threshold models
% sse is sum of squared errors, yfitted are the modelled values
[sse_thresh] = FitBrokenLine(b,x,y);

% get the slope for the H0 hypothesis that threshold = 0;
slp_linear = x\y;
sse_linear = sum((y - x.*slp_linear).^2);

% likelihood ratio
lr_statistic = length(x).*log(sse_linear/sse_thresh);

% degree of freedom lost in simpler model 
dof = 1;

% find the p value to tell whether the threshold model is indicated over the
% simpler linear model
p = chi2cdf(lr_statistic,dof);

% set outputs
thresh = b(1);
slope = b(2);
slope_linear = slp_linear;
% p_value<0.05 implies that the threshold exists
p_value = 1-p;

end

function [sse] = FitBrokenLine(b,x,y)
% b is a vector (x value of the threshold, slope of line)
% x and y are the arrays of point data

% unpack the variables to be fitted
x_thresh = b(1);
slope = b(2);

% get the estimated y values according to the fitted model
yfitted = zeros(size(y));
yfitted(x<x_thresh) = 0;
yfitted(x>=x_thresh) = slope.*(x(x>=x_thresh)-x_thresh);

% calculate the sum of squared errors between model and data
sse = sum((yfitted-y).^2);

end

%%%%%%%%%%%%%%%%%%%%%
%sub-function
%%%%%%%%%%%%%%%%%%%%%
function [f,xl,yl,tstat,yfitted] = FitBrokenLine2(b,x,y)
% if b defines the breakpoint of the x-values, fit a sloping line to y-data 
% below b, and fit a sloping line to the y-data with x above b.
% f is sum of squared errors, 
% (xl,yl) are points defining fitted line
% tstat is value of tstatistic for slope: b/(stderror of b)
% yfitted is value of fitted line at each x

f=0; % sum of squared errors between y-data and fitted broken line
bx = b(1);
by = b(2);
k=find(x<bx); n1=numel(k); % points below b
if ~isempty(k)
    % best fit line which must pass through (bx,by)
    q(1)=sum(y(k)-by)/sum(x(k)-bx); % q(1) is slope of best fit line
    q(2)=by-q(1)*bx; % q(2) is intercept of best fit line
    ndf=numel(x)-1; % ndf is number of points less than break point -1
    
    y2a=q(2)+q(1)*x(k); % calculate value of sloping line at each x below b
    SSEresid=sum((y(k)-y2a).^2)/(ndf-2); % sqrt of numerator of SE of slope coeff
    SXX=sum(x(k).^2)-(sum(x(k)).^2)/numel(x(k));% sqrt of denominator of SE of slope coeff
    Sb=sqrt(SSEresid)/sqrt(SXX); % Sb is SE of the slope coeff
    tstat=q(1)/Sb; % tscore of t-test that slope = 0
    f=f+sum((y(k)-y2a).^2); % sum of squared errors
    xl=[min(x) bx];yl=[q(2)+q(1)*[min(x) bx]]; % points at either end of sloping line
   yfitted(k)=y2a;

%     y1=mean(y(k)); %best fit flat line is the mean of y data
%     f=f+sum((y(k)-y1).^2); %sum of squared errors
%     xl=[min(x) b];yl=[y1 y1]; %points at either end of the fitted line
else % no data points to the left of threshold
    y1=[];
    xl=[];yl=[];
end

k=find(x>=bx); % points above b
if ~isempty(k)
    % find best fit line through rest of data
    if isempty(y2a) %if there was nothing to the left of threshold then we fit least squares
        p=polyfit(x,y,1);ndf=numel(x);
    else
        % best fit line which must pass through (b,y1)
        p(1)=sum(y(k)-by)/sum(x(k)-bx);p(2)=by-p(1)*bx;ndf=numel(x)-1;
    end
    y2=p(2)+p(1)*x(k); % calculate value of sloping line at each x above b
    SSEresid=sum((y(k)-y2).^2)/(ndf-2);
    SXX=sum(x(k).^2)-(sum(x(k)).^2)/numel(x(k));
    Sb=sqrt(SSEresid)/sqrt(SXX);tstat=p(1)/Sb;
    f=f+sum((y(k)-y2).^2); % sum of squared errors
    xl=[xl bx max(x)];yl=[yl p(2)+p(1)*[bx max(x)]]; % points at either end of sloping line
else
    % nothing to the right of threshold
    y2=[];tstat=[];
end
yfitted(k)=y2;

% statistic returned is a test of whether the slopes are equal

end

