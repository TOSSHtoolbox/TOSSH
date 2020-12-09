function [] = makeScatterPlot(CAMELS_signatures,CAMELS_data,lower_lim)
%makeScatterPlot Makes scatter plot to compare all CAMELS signatures.
%
%   INPUT
%   CAMELS_signatures: signatures calculated using TOSSH
%   CAMELS_data: CAMELS data
%   lower_lim: lower limit of color scale
%
%   OUTPUT
%   plot
%
%   Copyright (C) 2020
%   This software is distributed under the GNU Public License Version 3.
%   See <https://www.gnu.org/licenses/gpl-3.0.en.html> for details.

if nargin < 3
    error('Not enough input arguments.')
end

figure('pos',[100 100 750 450])

subplot(3,5,1); hold on; 
scatter(CAMELS_signatures.Q_mean,CAMELS_data.q_mean,25,CAMELS_data.flow_perc_complete,'filled')%,'k')
xlabel('Qmean TOSSH'); ylabel('Qmean CAMELS')
title(sprintf('\\rho_s = %.2f',corr(CAMELS_signatures.Q_mean,CAMELS_data.q_mean,'Type','Spearman','rows','complete')))
limits = [floor(min([min(CAMELS_signatures.Q_mean) min(CAMELS_data.q_mean)])) ...
    ceil(max([max(CAMELS_signatures.Q_mean) max(CAMELS_data.q_mean)]))]; 
xlim(limits); ylim(limits);
caxis([lower_lim 100]); %colorbar; 
plot(limits(1):.1:limits(2),limits(1):.1:limits(2),'k--');
% update cursor
dcm_obj = datacursormode(figure(gcf));
set(dcm_obj,'UpdateFcn',{@myupdatefcn,CAMELS_data.gauge_id,1:length(CAMELS_data.gauge_id)})
set(gca,'fontsize',7)

subplot(3,5,2); hold on; 
scatter(CAMELS_signatures.TotalRR,CAMELS_data.runoff_ratio,25,CAMELS_data.flow_perc_complete,'filled')%,'k')
xlabel('Q/P TOSSH'); ylabel('Q/P CAMELS')
title(sprintf('\\rho_s = %.2f',corr(CAMELS_signatures.TotalRR,CAMELS_data.runoff_ratio,'Type','Spearman','rows','complete')))
limits = [floor(min([min(CAMELS_signatures.TotalRR) min(CAMELS_data.runoff_ratio)])) ...
    ceil(max([max(CAMELS_signatures.TotalRR) max(CAMELS_data.runoff_ratio)]))]; 
xlim(limits); ylim(limits);
caxis([lower_lim 100]); %colorbar; 
plot(limits(1):.1:limits(2),limits(1):.1:limits(2),'k--');
set(gca,'fontsize',7)

subplot(3,5,3); hold on; 
scatter(CAMELS_signatures.QP_elasticity,CAMELS_data.stream_elas,25,CAMELS_data.flow_perc_complete,'filled')%,'k')
xlabel('QP elasticity TOSSH'); ylabel('QP elasticity CAMELS')
title(sprintf('\\rho_s = %.2f',corr(CAMELS_signatures.QP_elasticity,CAMELS_data.stream_elas,'Type','Spearman','rows','complete')))
limits = [floor(min([min(CAMELS_signatures.QP_elasticity) min(CAMELS_data.stream_elas)])) ...
    ceil(max([max(CAMELS_signatures.QP_elasticity) max(CAMELS_data.stream_elas)]))]; 
xlim(limits); ylim(limits);
caxis([lower_lim 100]); %colorbar; 
plot(limits(1):.1:limits(2),limits(1):.1:limits(2),'k--');
set(gca,'fontsize',7)
 
subplot(3,5,4); hold on; 
scatter(-CAMELS_signatures.FDC_slope,CAMELS_data.slope_fdc,25,CAMELS_data.flow_perc_complete,'filled')%,'k')
xlabel('FDC slope TOSSH'); ylabel('FDC slope CAMELS')
title(sprintf('\\rho_s = %.2f',corr(-CAMELS_signatures.FDC_slope,CAMELS_data.slope_fdc,'Type','Spearman','rows','complete')))
limits = [floor(min([min(-CAMELS_signatures.FDC_slope) min(CAMELS_data.slope_fdc)])) ...
    ceil(max([max(-CAMELS_signatures.FDC_slope) max(CAMELS_data.slope_fdc)]))]; 
xlim(limits); ylim(limits);
caxis([lower_lim 100]); %colorbar; 
plot(limits(1):.1:limits(2),limits(1):.1:limits(2),'k--');
set(gca,'fontsize',7)

subplot(3,5,5); hold on; 
scatter(CAMELS_signatures.BFI,CAMELS_data.baseflow_index,25,CAMELS_data.flow_perc_complete,'filled')%,'k')
xlabel('BFI TOSSH'); ylabel('BFI CAMELS')
title(sprintf('\\rho_s = %.2f',corr(CAMELS_signatures.BFI,CAMELS_data.baseflow_index,'Type','Spearman','rows','complete')))
limits = [floor(min([min(CAMELS_signatures.BFI) min(CAMELS_data.baseflow_index)])) ...
    ceil(max([max(CAMELS_signatures.BFI) max(CAMELS_data.baseflow_index)]))]; 
xlim(limits); ylim(limits);
caxis([lower_lim 100]); %colorbar; 
plot(limits(1):.1:limits(2),limits(1):.1:limits(2),'k--');
set(gca,'fontsize',7)

subplot(3,5,6); hold on; 
scatter(CAMELS_signatures.HFD_mean,CAMELS_data.hfd_mean,25,CAMELS_data.flow_perc_complete,'filled')%,'k')
xlabel('HFD TOSSH'); ylabel('HFD CAMELS')
title(sprintf('\\rho_s = %.2f',corr(CAMELS_signatures.HFD_mean,CAMELS_data.hfd_mean,'Type','Spearman','rows','complete')))
limits = [floor(min([min(CAMELS_signatures.HFD_mean) min(CAMELS_data.hfd_mean)])) ...
    ceil(max([max(CAMELS_signatures.HFD_mean) max(CAMELS_data.hfd_mean)]))]; 
xlim(limits); ylim(limits);
caxis([lower_lim 100]); %colorbar; 
plot(limits(1):.1:limits(2),limits(1):.1:limits(2),'k--');
set(gca,'fontsize',7)

subplot(3,5,7); hold on; 
scatter(CAMELS_signatures.Q5,CAMELS_data.q5,25,CAMELS_data.flow_perc_complete,'filled')%,'k')
xlabel('Q5 TOSSH'); ylabel('Q5 CAMELS')
title(sprintf('\\rho_s = %.2f',corr(CAMELS_signatures.Q5,CAMELS_data.q5,'Type','Spearman','rows','complete')))
limits = [floor(min([min(CAMELS_signatures.Q5) min(CAMELS_data.q5)])) ...
    ceil(max([max(CAMELS_signatures.Q5) max(CAMELS_data.q5)]))]; 
xlim(limits); ylim(limits);
caxis([lower_lim 100]); %colorbar; 
plot(limits(1):.1:limits(2),limits(1):.1:limits(2),'k--');
set(gca,'fontsize',7)

subplot(3,5,8); 
scatter(CAMELS_signatures.Q95,CAMELS_data.q95,25,CAMELS_data.flow_perc_complete,'filled')%,'k')
hold on;  
xlabel('Q95 TOSSH'); ylabel('Q95 CAMELS')
title(sprintf('\\rho_s = %.2f',corr(CAMELS_signatures.Q95,CAMELS_data.q95,'Type','Spearman','rows','complete')))
limits = [floor(min([min(CAMELS_signatures.Q95) min(CAMELS_data.q95)])) ...
    ceil(max([max(CAMELS_signatures.Q95) max(CAMELS_data.q95)]))]; 
xlim(limits); ylim(limits);
caxis([lower_lim 100]); %colorbar; 
plot(limits(1):.1:limits(2),limits(1):.1:limits(2),'k--');
set(gca,'fontsize',7)
 
subplot(3,5,9); hold on; 
scatter(CAMELS_signatures.high_Q_freq.*365,CAMELS_data.high_q_freq,25,CAMELS_data.flow_perc_complete,'filled')%,'k')
xlabel('High Q freq TOSSH'); ylabel('High Q freq CAMELS')
title(sprintf('\\rho_s = %.2f',corr(CAMELS_signatures.high_Q_freq.*365,CAMELS_data.high_q_freq,'Type','Spearman','rows','complete')))
limits = [floor(min([min(CAMELS_signatures.high_Q_freq.*365) min(CAMELS_data.high_q_freq)])) ...
    ceil(max([max(CAMELS_signatures.high_Q_freq.*365) max(CAMELS_data.high_q_freq)]))]; 
xlim(limits); ylim(limits);
caxis([lower_lim 100]); %colorbar; 
plot(limits(1):.1:limits(2),limits(1):.1:limits(2),'k--');
set(gca,'fontsize',7)

subplot(3,5,10); 
scatter(CAMELS_signatures.high_Q_dur,CAMELS_data.high_q_dur,25,CAMELS_data.flow_perc_complete,'filled')%,'k')
hold on; 
xlabel('High Q dur TOSSH'); ylabel('High Q dur CAMELS')
title(sprintf('\\rho_s = %.2f',corr(CAMELS_signatures.high_Q_dur,CAMELS_data.high_q_dur,'Type','Spearman','rows','complete')))
limits = [floor(min([min(CAMELS_signatures.high_Q_dur) min(CAMELS_data.high_q_dur)])) ...
    ceil(max([max(CAMELS_signatures.high_Q_dur) max(CAMELS_data.high_q_dur)]))]; 
xlim(limits); ylim(limits);
caxis([lower_lim 100]); %colorbar; 
plot(limits(1):.1:limits(2),limits(1):.1:limits(2),'k--');
set(gca,'fontsize',7)

subplot(3,5,11); hold on; 
scatter(CAMELS_signatures.low_Q_freq.*365,CAMELS_data.low_q_freq,25,CAMELS_data.flow_perc_complete,'filled')%,'k')
xlabel('Low Q freq TOSSH'); ylabel('Low Q freq CAMELS')
title(sprintf('\\rho_s = %.2f',corr(CAMELS_signatures.low_Q_freq.*365,CAMELS_data.low_q_freq,'Type','Spearman','rows','complete')))
limits = [floor(min([min(CAMELS_signatures.low_Q_freq.*365) min(CAMELS_data.low_q_freq)])) ...
    ceil(max([max(CAMELS_signatures.low_Q_freq.*365) max(CAMELS_data.low_q_freq)]))]; 
xlim(limits); ylim(limits);
caxis([lower_lim 100]); %colorbar; 
plot(limits(1):.1:limits(2),limits(1):.1:limits(2),'k--');
set(gca,'fontsize',7)

subplot(3,5,12); hold on; 
scatter(CAMELS_signatures.low_Q_dur,CAMELS_data.low_q_dur,25,CAMELS_data.flow_perc_complete,'filled')%,'k')
xlabel('Low Q dur TOSSH'); ylabel('Low Q dur CAMELS')
title(sprintf('\\rho_s = %.2f',corr(CAMELS_signatures.low_Q_dur,CAMELS_data.low_q_dur,'Type','Spearman','rows','complete')))
limits = [floor(min([min(CAMELS_signatures.low_Q_dur) min(CAMELS_data.low_q_dur)])) ...
    ceil(max([max(CAMELS_signatures.low_Q_dur) max(CAMELS_data.low_q_dur)]))]; 
xlim(limits); ylim(limits);
caxis([lower_lim 100]); %colorbar; 
plot(limits(1):.1:limits(2),limits(1):.1:limits(2),'k--');
set(gca,'fontsize',7)

subplot(3,5,13); hold on; 
scatter(CAMELS_signatures.zero_Q_freq,CAMELS_data.zero_q_freq,25,CAMELS_data.flow_perc_complete,'filled')%,'k')
xlabel('Zero Q freq TOSSH'); ylabel('Zero Q freq CAMELS')
title(sprintf('\\rho_s = %.2f',corr(CAMELS_signatures.zero_Q_freq,CAMELS_data.zero_q_freq,'Type','Spearman','rows','complete')))
limits = [floor(min([min(CAMELS_signatures.zero_Q_freq) min(CAMELS_data.zero_q_freq)])) ...
    ceil(max([max(CAMELS_signatures.zero_Q_freq) max(CAMELS_data.zero_q_freq)]))]; 
xlim(limits); ylim(limits);
caxis([lower_lim 100]); %colorbar; 
plot(limits(1):.1:limits(2),limits(1):.1:limits(2),'k--');
set(gca,'fontsize',7)

c = colorbar;
cmap = colormap;
colormap(flipud(cmap));
title(c,'% complete')
% x1=get(gca,'position');
% x=[0.63 0.11 0.01 0.20];
% x=[0.37 0.13 0.01 0.11];
x=[0.64 0.135 0.01 0.14];
set(c,'Position',x)
% set(gca,'position',x1)

% figure; hold on; 
% colormap(colour_mat);
% scatter(round(100.*results.P_mean)./100,CAMELS_data.p_mean,25,CAMELS_data.flow_perc_complete,'filled')%,'k')
% xlabel('Pmean TOSSH'); ylabel('Pmean CAMELS')
% title(sprintf('\\rho_s = %.2f',corr(results.P_mean,CAMELS_data.p_mean,'Type','Spearman','rows','complete')))
% limits = [floor(min([min(results.P_mean) min(CAMELS_data.p_mean)])) ...
%     ceil(max([max(results.P_mean) max(CAMELS_data.p_mean)]))]; 
% xlim(limits); ylim(limits);
% caxis([lower_lim 100]); %colorbar; 
% plot(limits(1):.1:limits(2),limits(1):.1:limits(2),'k--');
% % update cursor
% dcm_obj = datacursormode(figure(gcf));
% set(dcm_obj,'UpdateFcn',{@myupdatefcn,CAMELS_data.gauge_id,1:length(CAMELS_data.gauge_id)})

end

