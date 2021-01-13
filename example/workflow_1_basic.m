%% TOSSH workflow 1 - basic workflow (also shown in online documentation)
%
%   This script shows the basic functionalities of TOSSH with some example 
%   data. You can go through this example step by step by evaluating each 
%   section (separated by %%) using Ctrl+Enter. 
%
%   The example data used in this workflow are taken from CAMELS-GB 
%   (Coxon et al., 2020), see README_example_data.txt for more information
%   on data sources.
%
%   Copyright (C) 2020
%   This software is distributed under the GNU Public License Version 3.
%   See <https://www.gnu.org/licenses/gpl-3.0.en.html> for details.

close all
% clear all
clc

%% Check if required MATLAB toolboxes are installed
% Required are:
%   - MATLAB (TOSSH was developed using Matlab R2020a)
%   - Statistics and Machine Learning Toolbox
%   - Optimization Toolbox
% We can check which toolboxes we have installed using "ver".
% ver

%% Add directories to path
% We navigate to the TOSSH directory and add it to the Matlab path. 
% If we already are in this directory, we can use the pwd command:
mydir = pwd;
% Alternatively, we can specify my_dir manually:
% mydir = 'D:/Sebastian/Documents/MATLAB/TOSSH';
cd(mydir)
addpath(genpath(mydir));

%% Input data
% Every signature requires a streamflow (Q) time series and a corresponding 
% date vector (t). The date vector should be in datetime format, but 
% datenum format is also accepted and internally converted. Here is an 
% example of Q and t vectors in the correct format:
Q = [1.14;1.07;2.39;2.37;1.59;1.25;1.35;1.16;1.27;1.14]; 
t = [datetime(1999,10,1,0,0,0):datetime(1999,10,10,0,0,0)]';

% Typically, users will have their own data which they want to analyse. 
% We provide an example file to get a more realistic time series.
% The example file also contains precipitation (P), potential
% evapotranspiration (PET), and temperature (T) data, which are required 
% for some signatures. The paths are relative and assume that we are in 
% the TOSSH directory.
path = './example/example_data/'; % specify path
data = load(strcat(path,'33029_daily.mat')); % load data
t = data.t;
Q = data.Q; % streamflow [mm/day]
P = data.P; % precipitation [mm/day]
% Note: PET and T example data are provided but not used here.
% PET = data.PET; % potential evapotranspiration [mm/day]
% T = data.T; % temperature [degC]

%% Plot data
% We can plot the data to get a first idea of the hydrograph.
figure('pos',[100 100 350 200])
plot(t,Q,'k-','linewidth',1.0)
xlabel('Date')
ylabel('Streamflow [mm/day]')

% More information on the input data can be found here:
% https://TOSSHtoolbox.github.io/TOSSH/p1_overview.html.

%% Calculate signatures
% Once the input data are loaded, we can calculate different signatures.
% We start by calculating the mean flow Q_mean.
Q_mean = sig_Q_mean(Q,t);
fprintf('Q_mean = %.2f \n',Q_mean)

% Some signatures can be calculated using different methods and/or
% parameter values. For example, there are different options to calculate 
% the baseflow index (BFI). The default method is the UKIH smoothed minima 
% method with a parameter of 5 days.
BFI_UKIH = sig_BFI(Q,t);
fprintf('BFI_UKIH = %.2f \n',BFI_UKIH)
% Alternatively, we can use the Lyne-Hollick filter with a filter parameter 
% of 0.925.
BFI_LH = sig_BFI(Q,t,'method','Lyne_Hollick');
fprintf('BFI_LH = %.2f \n',BFI_LH)
% We can also change the parameter value of the UKIH method to 10 days.
BFI_UKIH10 = sig_BFI(Q,t,'method','UKIH','parameters',10);
fprintf('BFI_UKIH10 = %.2f \n',BFI_UKIH10)
% As we can see, all three options lead to slightly different values.
% More details and examples on the different methods/parameters can be
% found in the code of each function (e.g. sig_BFI.m).

% Some signatures also require precipitation (P) input time series. For
% example, the total runoff ratio requires both Q and P time series.
TotalRR = sig_TotalRR(Q,t,P);
fprintf('TotalRR = %.2f \n',TotalRR)

% Some signature functions come with a plotting functionality. For example, 
% we can calculate the slope of the flow duration curve (FDC) and plot the
% result.
[FDC_slope] = sig_FDC_slope(Q,t,'plot_results',true);
fprintf('FDC_slope = %.2f \n',FDC_slope)

% Some signatures are combinations of existing signatures, e.g. the
% baseflow fraction (K_b) defined as the ratio between mean baseflow Q_b  
% and mean precipitation P. This signature can also be calculated as
% K_b = BFI*TotalRR. We therefore do not provide an extra signature
% function, but suggest to use the two existing functions.
K_b = sig_BFI(Q,t)*sig_TotalRR(Q,t,P);
fprintf('K_b = %.2f \n',K_b)

% More information on the signatures contained in TOSSH can be found here:
% https://sebastiangnann.github.io/TOSSH_development/p2_signatures.html

%% Warnings and errors
% Each signature function can return a warning/error output. These 
% warning/error outputs indicate problems during signature calculation, but
% they do not stop code execution like a normal Matlab error would do. 
% Two outputs can be retrieved: an error flag (error_flag) that corresponds 
% to a certain type of warning/error, and a string (error_str) that 
% decribes the warning/error. If multiple warnings/errors occur, they are
% all listed in the error string, starting with the one that occurred last.

% A warning (error_flag = 1) typically indicates that the signature can be 
% calculated but should be interpreted with care, e.g. because there are 
% NaN values in the time series.
Q(1:10) = NaN;
[Q_mean, error_flag, error_str] = sig_Q_mean(Q,t);
fprintf('Q_mean = %.2f \n',Q_mean)
fprintf(error_str+"\n")
% We get the same mean value as before since the ten removed values do not
% influence the result much. In other cases, significant amounts of NaN
% entries might cause more problems.

% An error (error_flag > 1) indicates that the signature could not be 
% calculated, e.g. because there is a problem with the input data. For 
% example, if the input time series contains negative and thus physically 
% impossible values, NaN is returned. 
Q(1:10) = -1.0;
[Q_mean, error_flag, error_str] = sig_Q_mean(Q,t);
fprintf('Q_mean = %.2f \n',Q_mean)
fprintf(error_str+"\n")
% Since these warnings/errors do not stop the execution of the code, we can
% run the signature code for many catchments without breaking, even if
% for some of the catchments the signature cannot be calculated.

% There are also "normal" errors which can happen if the input parameters
% are specified incorrectly (wrong format, wrong range, etc.).
% For example, if we swap Q and t in the input, we will get an error. 
Q_mean = sig_Q_mean(t,Q);
% While this will stop code execution, such errors should be easily 
% avoidable by specifying all inputs correctly. 

%% Further information 
% Further information can be found in the online documentation: 
% https://TOSSHtoolbox.github.io/TOSSH/ and in the other example scripts
% (e.g. 'workflow_2_advanced.m').
