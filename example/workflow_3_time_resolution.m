%% TOSSH workflow 3 - comparison of time series with different resolution
%
%   This script shows how to use TOSSH with example data from the same 
%   catchment but with different time resolution. 
%
%   The example datasets used in this workflow were provided by the 
%   Environment Agency, see README_example_data.txt for more information on 
%   data sources.
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

% Daily, hourly, and 15min data from the same catchment.
data = load(strcat(path,'33029_multiple.mat')); % load data
t_daily = data.t_daily;
Q_daily = data.Q_daily; % streamflow [mm/day]
t_hourly = data.t_hourly;
Q_hourly = data.Q_hourly; % streamflow [mm/hour]
t_15min = data.t_15min;
Q_15min = data.Q_15min; % streamflow [mm/d15minay]

clear data

%% Plot data
% When we compare different time resolutions, we can see that some of the
% small peaks are smoothed out in the daily series.
figure('pos',[100 100 350 200])
hold on
plot(t_daily,Q_daily,'b')
plot(t_hourly,Q_hourly*24,'r-')
plot(t_15min,Q_15min*4*24,'g-.')
xlabel('Date')
ylabel('Streamflow [mm/day]')
xlim([datetime(2015,10,01) datetime(2016,09,30)])
legend('Daily','Hourly','15min','location','best')

%% Compare time resolution
% Different time resolutions have an impact on the calculation of
% hydrological signatures. We again calculate the BFI and the slope of the
% FDC, but now for the same catchment with different time resolutions.
BFI_daily = sig_BFI(Q_daily,t_daily);
BFI_hourly = sig_BFI(Q_hourly,t_hourly);
BFI_15min = sig_BFI(Q_15min,t_15min);
FDC_slope_daily = sig_FDC_slope(Q_daily,t_daily);
FDC_slope_hourly = sig_FDC_slope(Q_hourly,t_hourly);
FDC_slope_15min = sig_FDC_slope(Q_15min,t_15min);
RLD_daily = sig_RisingLimbDensity(Q_daily,t_daily,'eps',0.001*median(Q_daily,'omitnan'));
RLD_hourly = sig_RisingLimbDensity(Q_hourly,t_hourly,'eps',0.001*median(Q_hourly,'omitnan'));
RLD_15min = sig_RisingLimbDensity(Q_15min,t_15min,'eps',0.001*median(Q_15min,'omitnan'));
% We can store the results in a table to print them in the command window.
VarNames = {'BFI','FDC_slope','RLD'};
RowNames = {'1d','1h','15min'};
T = table([BFI_daily; BFI_hourly; BFI_15min],...
    [FDC_slope_daily; FDC_slope_hourly; FDC_slope_15min],...
    [RLD_daily; RLD_hourly*24; RLD_15min*(24*4)],...
    'VariableNames',VarNames,'RowNames',RowNames);
disp(T)
% From the results displayed in the command window and from the plots we
% can see that the BFI is sensitive to the timestep. This is because we
% always used the default method, which is the Lyne-Hollick filter with a
% parameter value of 0.925. The slope of the FDC is insensitive to the 
% timestep. No parameters need to be specified. 

% We also need to be careful with the units. For example, if we want to 
% compare mean flows we need to adjust the time series so that they have
% the same units (e.g. mm/day).
Q_mean_daily = sig_Q_mean(Q_daily,t_daily);
Q_mean_hourly = sig_Q_mean(Q_hourly,t_hourly);
Q_mean_15min = sig_Q_mean(Q_15min,t_15min);
VarNames = {'Q_mean','Q_mean_adj'};
T = table([Q_mean_daily; Q_mean_hourly; Q_mean_15min],...
    [Q_mean_daily; Q_mean_hourly*24; Q_mean_15min*4*24],...
    'VariableNames',VarNames,'RowNames',RowNames);
disp(T)

%% Further information 
% Further information can be found in the online documentation: 
% https://TOSSHtoolbox.github.io/TOSSH/ and in the other example scripts.
