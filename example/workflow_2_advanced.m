%% TOSSH workflow 2 - advanced workflow
%
%   This script shows how to use TOSSH with example data from multiple 
%   catchments. We first show how TOSSH can be used to analyse catchment
%   data. We then show typical issues that need to be taken care of, e.g. 
%   signatures that cannot be calculated for certain input data. 
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

%% Add directories to path
% We navigate to the TOSSH directory and add it to the Matlab path. 
% If we already are in this directory, we can use the pwd command:
mydir = pwd;
% Alternatively, we can specify my_dir manually:
% mydir = 'D:/Sebastian/Documents/MATLAB/TOSSH';
cd(mydir)
addpath(genpath(mydir));

%% Load data
path = './example/example_data/'; % specify path

% Daily data from 3 different catchments.
% Catchment 1
data = load(strcat(path,'33029_daily.mat')); % load data
t_1 = data.t;
Q_1 = data.Q; % streamflow [mm/day]
P_1 = data.P; % precipitation [mm/day]
% Note: PET and T example data are provided but not used here.
% PET_1 = data.PET; % potential evapotranspiration [mm/day]
% T_1 = data.T; % temperature [degC]
% Catchment 2
data = load(strcat(path,'39020_daily.mat')); % load data
t_2 = data.t;
Q_2 = data.Q; % streamflow [mm/day]
P_2 = data.P; % precipitation [mm/day]
% Catchment 3
data = load(strcat(path,'73014_daily.mat')); % load data
t_3 = data.t;
Q_3 = data.Q; % streamflow [mm/day]
P_3 = data.P; % precipitation [mm/day]

clear data

%% Plot data
% We can plot the data to get a first idea of the hydrographs. 
% The three different catchments show a different response, ranging from
% flashy to rather damped.
figure('pos',[100 100 350 200])
hold on
plot(t_1,Q_1./mean(Q_1),'b')
plot(t_2,Q_2./mean(Q_2),'r')
plot(t_3,Q_3./mean(Q_3),'g')
xlabel('Date')
ylabel('Q/mean(Q) [-]')
xlim([datetime(2002,10,01) datetime(2005,09,30)])
legend('Catchment 1','Catchment 2','Catchment 3','location','best')

%% Calculate signatures
% We can calculate some hydrological signatures for the three different
% catchments, such as the BFI or the slope of the flow duration.
BFI_1 = sig_BFI(Q_1,t_1,'method','UKIH');
BFI_2 = sig_BFI(Q_2,t_2,'method','UKIH');
BFI_3 = sig_BFI(Q_3,t_3,'method','UKIH');
FDC_slope_1 = sig_FDC_slope(Q_1,t_1);
FDC_slope_2 = sig_FDC_slope(Q_2,t_2);
FDC_slope_3 = sig_FDC_slope(Q_3,t_3);
% We can store the results in a table to print them in the command window.
VarNames = {'BFI','FDC_slope'};
RowNames = {'Catchment 1';'Catchment 2';'Catchment 3'};
Tab = table([BFI_1; BFI_2; BFI_3],...
    [FDC_slope_1; FDC_slope_2; FDC_slope_3],...
    'VariableNames',VarNames,'RowNames',RowNames);
disp(Tab)
% The signatures (shown in the printed table) correspond well with the 
% hydrographs. Catchments 1 and 2 show a damped response indicated by high
% BFIs and low slopes of the FDC, while catchment 3 shows a flashy response
% indicated by a lower BFI and a higher slope of the FDC.

%% Calculation functions
% We can calculate all signatures from one of the signature sets using a
% calculation function, e.g. the Euser et al. (2013) calculation fucntion.
% First, we need to create cell arrays containing the time series. 
Q_cell = {Q_1; Q_2; Q_3};
t_cell = {t_1; t_2; t_3};
Euser_signatures = calc_Euser(Q_cell,t_cell);
% The signatures are returned as a single struct file. We can then extract
% individual signatures for all three catchments, for example the lag-1
% autocorrelation AC1.
VarNames = {'AC1'};
Tab = table(Euser_signatures.AC1,...
    'VariableNames',VarNames,'RowNames',RowNames);
disp(Tab)

% Cell arrays are used to be able to handle time series of different 
% length which can be the case with real data sets. 
Q_cell = {Q_1(1:1000); Q_2(1:500); Q_3(2000:end)};
t_cell = {t_1(1:1000); t_2(1:500); t_3(2000:end)};
Euser_signatures_short = calc_Euser(Q_cell,t_cell);
% The shorter time series result in slightly different signature values.
VarNames = {'AC1'};
Tab = table(Euser_signatures_short.AC1,...
    'VariableNames',VarNames,'RowNames',RowNames);
disp(Tab)

%% Signatures with plotting functionality
% Some signatures come with a plotting functionality. This can help to
% determine the suitability of input parameters or to determine the 
% suitability of signature assumptions.
% If we perform a recession analysis and choose to plot the result, we 
% obtain a plot that shows the recession segments and a plot that shows the 
% the fitted recession curves coloured according to the season. 
[para_mat, ~, ~, ~, fig_RecessionAnalysis] = ...
    sig_RecessionAnalysis(Q_3,t_3,'plot_results',true);
% We can also obtain the figure handles which allows us to manipulate that
% figure. For example, we can zoom into the time series to inspect the
% chosen recession segments.
figure(fig_RecessionAnalysis.RecessionSegments)
xlim([datetime(2002,10,01) datetime(2003,03,31)])

% We can also calculate the recession constant using a master recession 
% curve which assumes exponential recession behaviour and check whether
% this assumption is reasonable. Note that we set the eps parameter to a
% small value that allows for recessions with small "bumps".
[K, ~, ~, fig_BaseflowRecessionK] = ...
    sig_BaseflowRecessionK(Q_3,t_3,'plot_results',true,'eps',0.01*median(Q_1));

% Some signatures focus on patterns that are best analysed visually. For 
% example, the quickflow response during a rainfall event can indicate 
% which runoff generation processes occur in a catchment. 
[ie_effect, se_effect, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, fig_EventGraphThresholds] = ...
    sig_EventGraphThresholds(Q_3,t_3,P_3,'plot_results',true);
% We can see that most summer storms lead to little runoff, while all 
% winter storms produce runoff. This indicates saturation excess quickflow.

%% Signatures that are variations of the same signature
% Some signatures are variations of existing signatures. For example, many
% papers use differen sections of the flow duration curve. The default 
% option is the slope between the 33th and the 66th percentile.
FDC_midslope = sig_FDC_slope(Q_1,t_1,'plot_results',true); 
% Alternatively we can calculate the high-section or the low-section slope
% (e.g. Yilmaz et al., 2008). Note that the percentiles are defined as 
% exceedance probabilities, i.e. the high-section corresponds to low values
% (e.g. [0 0.3]) and vice versa.
FDC_highslope = sig_FDC_slope(Q_1,t_1,'slope_range', [0 0.3], 'plot_results',true); 
FDC_lowslope = sig_FDC_slope(Q_1,t_1,'slope_range', [0.8 1.0], 'plot_results',true);

%% NaN values
% A common problem are NaN values (or missing values) in time series. TOSSH
% can handle NaN values, for example by calculating the signature using the
% remaining series. Many missing values can influence the resulting 
% signature value. To show the effect of NaN values we can create a new
% flow time series that contains NaN values.
Q_mean_1 = sig_Q_mean(Q_1,t_1);
Q_tmp = Q_1;
Q_tmp(1:50) = NaN;
Q_mean_2 = sig_Q_mean(Q_tmp,t_1);
Q_tmp(1:2000) = NaN;
Q_mean_3 = sig_Q_mean(Q_tmp,t_1);
VarNames = {'Q_mean'};
RowNames = {'Full series','A few NaNs','Many NaNs'};
Tab = table([Q_mean_1; Q_mean_2; Q_mean_3],...
    'VariableNames',VarNames,'RowNames',RowNames);
disp(Tab)
% We can see that few missing values do not change the mean flow, but many
% NaN values have an effect on the signature. 

%% Utility functions
% Utility functions are used internally but are coded separately to enhance
% reusability, both within TOSSH and externally. For example, we can use
% our recession extraction utility function to extract recession segments.
flow_section = util_RecessionSegments(Q_1(1:365), t_1(1:365), ...
    'recession_length',5,'start_of_recession','peak','plot_results',true);
% util_RecessionSegments returns indices corresponding to the start and end
% of each recession segment.

%% Warnings and errors
% NaN values or other problems during signature calculation can be detected 
% by looking at the warning/error outputs each signature functions returns. 
% Two outputs can be retrieved: an error flag (error_flag) that corresponds  
% to a certain type of warning/error, and a string (error_str) that 
% decribes the warning/error (see also workflow_1_basic.m). If 
% multiple warnings/errors occur, they are all listed in the error string, 
% starting with the one that occurred last.

% Some signatures cannot be calculated for certain catchments. For example,
% if a catchment is highly ephemeral, the slope of the flow duration curve
% cannot be calculated. We can simulate such a catchment by setting many
% values to 0.
Q_tmp = Q_1;
Q_tmp(1:2000) = 0;
[FDC_slope, error_flag, error_str] = sig_FDC_slope(Q_tmp,t_1);
fprintf('FDC_slope = %.2f \n',FDC_slope)
fprintf(error_str+"\n")
% The function returns NaN and the error statement indicates that the slope 
% of the FDC could not be calculated.

% There are also "normal" errors which can happen if the input parameters
% are specified incorrectly (wrong format, wrong range, etc.).
[FDC_slope, error_flag, error_str] = sig_FDC_slope(Q_tmp,t_1,'slope_range',[0.33 1.1]);

% Note that we can turn off all Matlab warnings using warning('off','all').

%% Further information 
% Further information can be found in the online documentation: 
% https://TOSSHtoolbox.github.io/TOSSH/ and in the other example scripts.
