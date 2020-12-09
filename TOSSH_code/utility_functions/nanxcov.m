function [xycov,lags,nanp] = nanxcov(x,y,option1,option2)

% NANXCOV   Cross-covariance function estimates. Skips NaNs.
%   XCOV(A,B), where A and B are length M vectors, returns the
%   length 2*M-1 cross-covariance sequence in a column vector.
%
%   XCOV(A), when A is a vector, is the auto-covariance sequence.
%   XCOV(A), when A is an M-by-N matrix, is a large matrix with
%   2*M-1 rows whose N^2 columns contain the cross-covariance
%   sequences for all combinations of the columns of A.
%   The zeroth lag of the output covariance is in the middle of the
%   sequence, at element or row M.
%
%   The cross-covariance is the cross-correlation function of
%   two sequences with their means removed:
%        C(m) = E[(A(n+m)-MA)*conj(B(n)-MB)]
%   where MA and MB are the means of A and B respectively.
%
%   XCOV(...,MAXLAG) computes the (auto/cross) covariance over the
%   range of lags:  -MAXLAG to MAXLAG, i.e., 2*MAXLAG+1 lags.
%   If missing, default is MAXLAG = M-1.
%
%   [C,LAGS] = XCOV(...) returns a vector of lag indices (LAGS).
%
%   XCOV(...,SCALEOPT), normalizes the covariance according to SCALEOPT:
%       biased   - scales the raw cross-covariance by 1/M.
%       unbiased - scales the raw covariance by 1/(M-abs(k)), where k
%                  is the index into the result.
%       coeff    - normalizes the sequence so that the covariances at
%                  zero lag are identically 1.0.
%       none     - no scaling (this is the default).
%
%   See also XCORR, CORRCOEF, CONV, COV and XCORR2.

%   Author(s): L. Shure, 1-9-88



% skips NaNs, P. Sturm, 10-Nov-2008
% outputs isnan vector of skipped NaNs, P. Sturm, 9-Nov-2009

%   References:
%     [1] J.S. Bendat and A.G. Piersol, "Random Data:
%         Analysis and Measurement Procedures", John Wiley
%         and Sons, 1971, p.332.
%     [2] A.V. Oppenheim and R.W. Schafer, Digital Signal
%         Processing, Prentice-Hall, 1975, pg 539.

nanx = isnan(x);
nany = isnan(y);
if nargin == 1
    x(nanx) = [];
    mx = size(x, 1);
    [xycov,l] = xcorr(x-ones(mx,1)*mean(x));
    nanp = nanx;
elseif nargin == 2
    x(nanx|nany) = [];
    y(nanx|nany) = [];
    mx = size(x, 1);
    if ischar(y)||(~ischar(y)&&length(y)==1)
        [xycov,l] = xcorr(x-ones(mx,1)*mean(x),y);
    else
        my = size(y, 1);
        [xycov,l] = xcorr(x-ones(mx,1)*mean(x),y-ones(my,1)*mean(y));
    end
    nanp = (nanx|nany);
elseif nargin == 3
    x(nanx|nany) = [];
    y(nanx|nany) = [];
    mx = size(x, 1);
    my = size(y, 1);
    if length(y)==1
        [xycov,l] = xcorr(x-ones(mx,1)*mean(x),y,option1);
    else
        [xycov,l] = xcorr(x-ones(mx,1)*mean(x),y-ones(my,1)*mean(y),option1);
    end
    nanp = (nanx|nany);
elseif nargin == 4
    x(nanx|nany) = [];
    y(nanx|nany) = [];
    mx = size(x, 1);
    my = size(y, 1);
    [xycov,l] = xcorr(x-ones(mx,1)*mean(x),y-ones(my,1)*mean(y),...
                         option1,option2);
    nanp = (nanx|nany);
end
if nargout > 1
    lags = l;
end