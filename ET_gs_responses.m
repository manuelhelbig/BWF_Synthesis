%% plotting ET and gs response to VPD %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% created by M. Helbig (helbigm@mcmaster.ca) on 2020-02-04
% used for Helbig et al. (????) Increasing contribution of peatlands to boreal evapotranspiration in a warming climate
function ET_gs_responses(PTL, ET_BINS_UP, GS_BINS_UP, VPD_BINS, par_gs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input:
% PTL: 0-1 vector flagging peatland (1) vs. forest (0) sites
% ET_BINS_UP: matrix of upper boundaries of evapotranspiration per VPD bin across all sites [mm hr-1]
% GS_BINS_UP: matrix of upper boundaries of surface conductance per VPD bin across all sites [m s-1]
% VPD_BINS: matrix of medians per VPD bin across all sites
% par_gs: model parameters for surface conductance response for all sites
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot gs and ET (upper boundary)
% calculate probability of ETPTL > ETFOR
i_PTL = find(PTL==1);
i_FOR = find(PTL==0);
% gs function
gsFxn=@(params,VPD) (params(2)+(1+params(1)./(sqrt(VPD))));
% VPD bins
PRC_VPD=0.1:0.1:3;
for p=1:length(PRC_VPD)+1;
    
    x_up=ET_BINS_UP(i_PTL,p);
    y_up=ET_BINS_UP(i_FOR,p);

    x1_up=x_up(~isnan(x_up));
    y1_up=y_up(~isnan(y_up));
    % get the permutationTest function here:
    % https://www.mathworks.com/matlabcentral/fileexchange/63276-permutation-test
    [pval_ET(p), observeddifference(p), effectsize] = permutationTest(x1_up,y1_up,999,'sidedness','larger');
    clear x y x1 y1 x_up y_up x1_up y1_up
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%

map_CB=flipud(brewermap([9],'BrBg'));
% ET response to VPD
figure,
subplot(2,2,1);

% plot peatland response
% allocate vector for standard error estimates
ET_SE_UP=NaN(1,length(PRC_VPD)+1);
for k=1:length(PRC_VPD)+1
    % standard error per bin (standard deviation/number of sites)
    ET_SE_UP(k)=nanstd(ET_BINS_UP(i_PTL,k))./sqrt(sum(~isnan(ET_BINS_UP(i_PTL,k))));
end
hold on
% use threshold below which data is not shown
thr = 0;
VPD_filt=find(nanmean(VPD_BINS(i_PTL,:))>thr);
% plot standard error of peatland response
plotshaded(nanmean(VPD_BINS(i_PTL,VPD_filt)),[(nanmean(ET_BINS_UP(i_PTL,VPD_filt))-ET_SE_UP(VPD_filt))' (nanmean(ET_BINS_UP(i_PTL,VPD_filt))+ET_SE_UP(VPD_filt))'],'r');
hold on
% plot mean peatland response
l3=plot(nanmean(VPD_BINS(i_PTL,VPD_filt)),nanmean(ET_BINS_UP(i_PTL,VPD_filt)),':','Color',map_CB(2,:),'LineWidth',2);

% plot forest response
% allocate vector for standard error estimates
ET_SE_UP=NaN(1,length(PRC_VPD)+1);
for k=1:length(PRC_VPD)+1
    % standard error per bin (standard deviation/number of sites)
    ET_SE_UP(k)=nanstd(ET_BINS_UP(i_FOR,k))./sqrt(sum(~isnan(ET_BINS_UP(i_FOR,k))));
end
plotshaded(nanmean(VPD_BINS(i_FOR,VPD_filt)),[(nanmean(ET_BINS_UP(i_FOR,VPD_filt))-ET_SE_UP(VPD_filt))' (nanmean(ET_BINS_UP(i_FOR,VPD_filt))+ET_SE_UP(VPD_filt))'],'r');
l4=plot(nanmean(VPD_BINS(i_FOR,VPD_filt)),nanmean(ET_BINS_UP(i_FOR,VPD_filt)),':','Color',map_CB(8,:),'LineWidth',2);

axis square
legend([l3 l4],{'Peatlands','Forests'});
% set x-axis limits
xlim([-0.1 3.8]);
% add axis labels
ylabel('ET [mm hr^{-1}]');
xlabel('VPD [kPa]');
set(gca,'FontSize',16);

% plot difference between peatland and forest ET
subplot(2,2,2);
for k = 1:1000;
    % randomly sample with replacement (bootstrap)
    i_PTL_shuffle=randsample(length(i_PTL),length(i_PTL),1);
    i_FOR_shuffle=randsample(length(i_FOR),length(i_FOR),1);
    DeltaET_UP(:,k)=nanmean(ET_BINS_UP(i_PTL(i_PTL_shuffle),:))-nanmean(ET_BINS_UP(i_FOR(i_FOR_shuffle),:));
end
% plot mean difference
plot(nanmean(VPD_BINS(i_PTL,VPD_filt)),nanmean(DeltaET_UP(VPD_filt,:)'),'-','LineWidth',2);
hold on
% color code according to permutation test results
scatter(nanmean(VPD_BINS(:,VPD_filt)),nanmean(DeltaET_UP(VPD_filt,:)'),35,(1-pval_ET(VPD_filt)),'filled');
cbr=colorbar;
caxis([0.9 1]);
set(cbr,'YTickLabel',(0.9:0.05:1).*100)
set(get(cbr,'ylabel'),'String', {'Probability of ET_{PTL} > ET_{FOR} [%]'},'FontSize',14);
colormap(brewermap([10],'GnBu'))

ylabel('\DeltaET [mm hr^{-1}]');
xlabel('VPD [kPa]');
set(gca,'FontSize',16);

hold on

axis square
xlim([-0.1 3.8]);

% plot surface conductance response
subplot(2,2,3)
sel=find(nanmean(VPD_BINS(i_PTL,:))>0.5);

GS_SE_UP=NaN(length(PRC_VPD)+1);
for k=1:length(PRC_VPD)+1
    GS_SE_UP(k)=(nanstd(GS_BINS_UP(i_PTL,k))./sqrt(sum(~isnan(GS_BINS_UP(i_PTL,k))))).*1000;
end
plotshaded(nanmean(VPD_BINS(i_PTL,sel)),[(nanmean(GS_BINS_UP(i_PTL,sel)).*1000-GS_SE_UP(sel))' (nanmean(GS_BINS_UP(i_PTL,sel)).*1000+GS_SE_UP(sel))'],'r');
hold on
l1=plot(nanmean(VPD_BINS(i_PTL,sel)),nanmean(GS_BINS_UP(i_PTL,sel)).*1000,'-','Color',map_CB(2,:),'LineWidth',3);

sel=find(nanmean(VPD_BINS(i_FOR,:))>0.5);

GS_SE_UP=NaN(length(PRC_VPD)+1);
for k=1:length(PRC_VPD)+1
    GS_SE_UP(k)=(nanstd(GS_BINS_UP(i_FOR,k))./sqrt(sum(~isnan(GS_BINS_UP(i_FOR,k))))).*1000;
end
plotshaded(nanmean(VPD_BINS(i_FOR,sel)),[(nanmean(GS_BINS_UP(i_FOR,sel)).*1000-GS_SE_UP(sel))' (nanmean(GS_BINS_UP(i_FOR,sel)).*1000+GS_SE_UP(sel))'],'r');
l2=plot(nanmean(VPD_BINS(i_FOR,sel)),nanmean(GS_BINS_UP(i_FOR,sel)).*1000,'-','Color',map_CB(8,:),'LineWidth',3);
axis square
% plot modelled surface conductance responses (using mean parameters for peatlands and forests)
l3=plot(nanmean(VPD_BINS(i_PTL,sel)),gsFxn(nanmean(par_gs(i_PTL,:)),nanmean(VPD_BINS(i_PTL,sel))),'--','Color',map_CB(3,:),'LineWidth',2);
hold on
l4=plot(nanmean(VPD_BINS(i_FOR,sel)),gsFxn(nanmean(par_gs(i_FOR,:)),nanmean(VPD_BINS(i_FOR,sel))),'--','Color',map_CB(7,:),'LineWidth',2);

xlim([0 3.5]);
legend([l1 l2 l3 l4],{'Peatlands','Forests','modelled PTL','modelled FOR'});
ylabel('g_s [mm s^{-1}]');
xlabel('VPD [kPa]');
set(gca,'FontSize',16);

% calculate if gs per VPD is significantly different
for p=1:length(PRC_VPD)+1;
    x=GS_BINS_UP(i_PTL,p);
    y=GS_BINS_UP(i_FOR,p);

    x1=x(~isnan(x));
    y1=y(~isnan(y));

    [pval_GS(p), observeddifference(p), effectsize] = permutationTest(x1,y1,999,'sidedness','larger');   
    clear x y x1 y1
end

for k = 1:1000;
    i_PTL_shuffle=randsample(length(i_PTL),length(i_PTL),1);
    i_FOR_shuffle=randsample(length(i_FOR),length(i_FOR),1);
    DeltaGS(:,k)=nanmean(GS_BINS_UP(i_PTL(i_PTL_shuffle),:))-nanmean(GS_BINS_UP(i_FOR(i_FOR_shuffle),:));
end
% plot difference in surface conductance
subplot(2,2,4)
sel=find(nanmean(VPD_BINS(i_PTL,:))>0.5);

plot(nanmean(VPD_BINS(i_PTL,sel)),nanmean(DeltaGS(sel,:)').*1000,'-','LineWidth',3);
hold on
scatter(nanmean(VPD_BINS(i_PTL,sel)),nanmean(DeltaGS(sel,:)').*1000,35,(1-pval_GS(sel)),'filled');
cbr=colorbar;
caxis([0.9 1]);
set(cbr,'YTickLabel',(0.9:0.05:1).*100)
set(get(cbr,'ylabel'),'String', {'Probability of g_{s {PTL}} > g_{s {FOR}}[%]'},'FontSize',14);
colormap(brewermap([10],'GnBu'))

ylabel('\Delta g_s [mm s^{-1}]');
xlabel('VPD [kPa]');
set(gca,'FontSize',16);
axis square
xlim([-0.1 3.5]);
