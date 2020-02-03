%% boundary line analysis %%
function [par] = boundary_line_analysis(lng, DOY, HH, H, LE, VPD, RN, SWIN, TA, GS)
%%% input variables %%%
% lng = number of sites
% DOY = day of year of dataset
% HH = timestamp (hours)
% H = sensible heat flux (W m-2)
% LE = latent heat flux (W m-2)
% VPD = vapour pressure deficit (kPa)
% RN = net radiation (W m-2)
% SWIN = incoming shortwave radiation (W m-2)
% TA = air temperature (degC)
% GS = surface conductance (m s-1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time_start = 15.5;
time_end = 18;
for k=1:lng;
    k
    % select surface conductance to be analysed
        % data between May and September
        % between 15.5 and 18h local time
        % minimum LE = 50 Wm-2
        % minimum available energy = 100 Wm-2
        
    sel = find((DOY<DOY_start | DOY>=DOY_end) | LE<50 ...
            | (HH<time_start | HH>=time_end) | (H+LE)<100);
    sel1 = find((DOY<DOY_start | DOY>=DOY_end));
    sel2 = find((DOY<DOY_start | DOY>=DOY_end) ...
            | (HH<time_start | HH>=time_end));
    
    % available energy (sensible + latent heat flux)
    AV=LE+H;
    
    ET = LE.*30.*60./Lv.*2;
    VPD(rem1)=NaN;
    GS(rem)=NaN;
    RN(rem)=NaN;
    AV(rem)=NaN;
    ET(rem2)=NaN;
    SWIN(rem)=NaN;
    TA(rem)=NaN;
    GS_FIT=GS;
    VPD_FIT=VPD;
    
    if sum(~isnan(GS))>10;
       GS_FIT(VPD_FIT<1)=NaN;
       VPD_FIT(VPD_FIT<1)=NaN;
       %x = log(VPD_FIT);
       x = VPD_FIT;
       y = GS_FIT;
       i = find(~isnan(x) & ~isnan(y));
       
       gsFxn=@(params,VPD) (params(1)+params(2).*log(VPD));
       par0=[-0.002 0.005];
       %par(k,:) = nlinfit(x,y,gsFxn,par0); %,opts);
       %par(k,:) = polyfit(x(i),y(i),1);
    else
       %par(k,:) = [NaN NaN];
       %par_gs(k,:) = [NaN];
    end   
    clear VPD_FIT GS_FIT x y
    % get aerodynamic conductance
    ga=VPD_LE{k}(:,10);
    ga(rem)=NaN;
    ga_SIT(k)=nanmedian(ga);
    
    % derive bins for VPD
    PRC_VPD=0.1:0.1:3;
    i=find(~isnan(GS));
    i_ET=find(~isnan(ET));
    VPD_UP = [];
    GS_UP = [];
    Rn_UP = [];
    for l =1:length(PRC_VPD)+1;
            if l == 1;
                ind=find(VPD<=PRC_VPD(l));     
            elseif l<=length(PRC_VPD);
                ind=find(VPD>PRC_VPD(l-1) & VPD<=PRC_VPD(l));  
            else
                ind=find(VPD>PRC_VPD(l-1));
            end
            
            if length(ind)<3 | isempty(i)
            	GS_BINS(k,l)=NaN;
                GS_BINS_SD(k,l)=NaN;  
                y1=NaN;
                
                VPD_sel = VPD(ind);
                RN_sel = RN(ind);
                AV_sel = AV(ind);
                ET_sel = ET(ind);

                % select outliers using quartile method
                IQR = prctile(ET_sel,75)-prctile(ET_sel,25);
                outl = find(ET_sel>prctile(ET_sel,75)+1.5*IQR | ET_sel<prctile(ET_sel,25)-1.5*IQR);
                ET_sel(outl)=NaN;
                
                VPD_BINS(k,l)=nanmean(VPD_sel);
                RN_BINS(k,l)=nanmean(RN_sel);
                ET_BINS(k,l)=nanmean(ET_sel);
                ET_BINS_SD(k,l)=nanstd(ET_sel);
                AV_BINS(k,l)=nanmean(AV_sel);
                VPD_BINS_SD(k,l)=nanstd(VPD_sel);
                
                % upper boundary
                x1=VPD_sel(ET_sel>=ET_BINS(k,l)+1.*ET_BINS_SD(k,l));
                z1=ET_sel(ET_sel>=ET_BINS(k,l)+1.*ET_BINS_SD(k,l));
                w1=RN_sel(ET_sel>=ET_BINS(k,l)+1.*ET_BINS_SD(k,l));
                v1=AV_sel(ET_sel>=ET_BINS(k,l)+1.*ET_BINS_SD(k,l));
                y1=NaN(size(x1));    
                clear GS_sel VPD_sel RN_sel
                VPD_BINS_UP(k,l)=nanmedian(x1);
                GS_BINS_UP(k,l)=nanmedian(y1);
                ET_BINS_UP(k,l)=nanmedian(z1);
                AV_BINS_UP(k,l)=nanmedian(v1);
                RN_BINS_UP(k,l)=nanmedian(w1);
                
            else
                VPD_sel = VPD(ind);
                RN_sel = RN(ind);
                GS_sel = GS(ind);
                AV_sel = AV(ind);
                ET_sel = ET(ind);
                % select outliers using quartile method
                IQR = prctile(GS_sel,75)-prctile(GS_sel,25);
                outl = find(GS_sel>prctile(GS_sel,75)+1.5*IQR | GS_sel<prctile(GS_sel,25)-1.5*IQR);
                GS_sel(outl)=NaN;
                
                % select outliers using quartile method
                IQR = prctile(ET_sel,75)-prctile(ET_sel,25);
                outl = find(ET_sel>prctile(ET_sel,75)+1.5*IQR | ET_sel<prctile(ET_sel,25)-1.5*IQR);
                %ET_sel(outl)=NaN;
%                 subplot(1,2,1)
%                 plot(VPD_sel,GS_sel,'.','Color',[0.6 0.6 0.6]);
%                 hold on
%                 axis square
%                 xlim([0 4])
%                 subplot(1,2,2)
%                 plot(VPD_sel,ET_sel,'.','Color',[0.6 0.6 0.6]);
%                 hold on
%                 xlim([0 4])
%                 axis square
                VPD_BINS(k,l)=nanmean(VPD_sel);
                RN_BINS(k,l)=nanmean(RN_sel);
                ET_BINS(k,l)=nanmean(ET_sel);
                GS_BINS(k,l)=nanmean(GS_sel);
                AV_BINS(k,l)=nanmean(AV_sel);
                VPD_BINS_SD(k,l)=nanstd(VPD_sel);
                GS_BINS_SD(k,l)=nanstd(GS_sel);
                ET_BINS_SD(k,l)=nanstd(ET_sel);
                % upper boundary
                x1=VPD_sel(GS_sel>=GS_BINS(k,l)+1.*GS_BINS_SD(k,l));
                y1=GS_sel(GS_sel>=GS_BINS(k,l)+1.*GS_BINS_SD(k,l));
                z1=ET_sel(GS_sel>=GS_BINS(k,l)+1.*GS_BINS_SD(k,l));
                w1=RN_sel(GS_sel>=GS_BINS(k,l)+1.*GS_BINS_SD(k,l));
                v1=AV_sel(GS_sel>=GS_BINS(k,l)+1.*GS_BINS_SD(k,l));
%                 subplot(1,2,1)
%                 plot(x1,y1,'o','Color',[0.1 0.1 0.1]);
%                 subplot(1,2,2)
%                 plot(x1,z1,'o','Color',[0.1 0.1 0.1]);
%                 x1=VPD_sel(ET_sel>=ET_BINS(k,l)+1.*ET_BINS_SD(k,l));
%                 y1=GS_sel(ET_sel>=ET_BINS(k,l)+1.*ET_BINS_SD(k,l));
%                 z1=ET_sel(ET_sel>=ET_BINS(k,l)+1.*ET_BINS_SD(k,l));
%                 w1=RN_sel(ET_sel>=ET_BINS(k,l)+1.*ET_BINS_SD(k,l));

                clear GS_sel VPD_sel RN_sel
                VPD_BINS_UP(k,l)=nanmedian(x1);
                GS_BINS_UP(k,l)=nanmedian(y1);
                ET_BINS_UP(k,l)=nanmedian(z1);
                AV_BINS_UP(k,l)=nanmedian(v1);
                RN_BINS_UP(k,l)=nanmedian(w1);
%                 subplot(1,2,1)
%                 plot(nanmedian(x1),nanmedian(y1),'o','Color',[0.1 0.6 0.1],'LineWidth',2);
%                 subplot(1,2,2)
%                 plot(nanmedian(x1),nanmedian(z1),'o','Color',[0.1 0.6 0.1],'LineWidth',2);
            end
            VPD_UP=vertcat(VPD_UP,x1);
            GS_UP=vertcat(GS_UP,y1);
            Rn_UP=vertcat(Rn_UP,w1);
            clear x1 y1 w1

    end 
    % fit log model for VPD-gs relationship
    
            if sum(~isnan(GS_UP))>10;
                % only fit model to VPD > 1kPa
                thr=1;
                GS_UP(VPD_UP<thr)=NaN;
                VPD_UP(VPD_UP<thr)=NaN;
                Rn_UP(VPD_UP<thr)=NaN;
                %x = log(VPD_UP);
                x1 = VPD_UP;
                %x1=horzcat(x1,Rn_UP);
                y = GS_UP.*1000;
                z=Rn_UP;
                %i = find(~isnan(x) & ~isnan(y));
                i1 = find(~isnan(x1) & ~isnan(y));
                %gsFxn=@(params,VPD) (params(1).*log(VPD));
                gsFxn=@(params,VPD) (params(2)+(1+params(1)./(sqrt(VPD))));
                
                %gsFxn=@(params,VPD) (params(3)+(1+params(1)./(VPD.^(params(2)))));
                %par0=[0.5 2];
                par0_gs=[10 -5];
                %par0_gs=[-5 0.5 10];
                par_gs(k,:) = nlinfit(x1(i1),y(i1),gsFxn,par0_gs,opts);
                %par(k,:) = polyfit(x(i),y(i),1);
                y=GS_BINS_UP(k,:);
                x=VPD_BINS_UP(k,:);
                RMSE(k)=sqrt(mean(y(~isnan(x) & ~isnan(y) & x>1)-(gsFxn(par_gs(k,:),x(~isnan(x) &~isnan(y) & x>1))./1000)).^2);
                R=corrcoef(y(~isnan(x) & ~isnan(y) & x>1),(gsFxn(par_gs(k,:),x(~isnan(x) &~isnan(y) & x>1))./1000));
                R2_VPD(k)=R(1,2)^2;
             
            else
                %par(k,:) = [NaN NaN];
                RMSE(k)=NaN;
                par_gs(k,:) = [NaN NaN];
                R2_VPD(k)=NaN;
            end    
            
            VPD_BOUND{k}=VPD_UP;
            GS_BOUND{k}=GS_UP;
            clear GS_UP VPD_UP

             if isempty(GS_BOUND{k}) | sum(~isnan(GS_BOUND{k}))==0
                 r=r+1;
                 NO_GS{r}=SITE_NAME{k};
             end
% derive Bowen ratio
    sel = find((DOY<DOY_start | DOY>=DOY_end) | VPD_LE{k}(:,2)<50 ...
            | (HH<time_start | HH>=time_end) | (VPD_LE{k}(:,8)+VPD_LE{k}(:,2))<100);
   
    VPD_GS=VPD_LE{k}(:,1);
    VPD_GS(sel)=NaN;
    Rnet_GS=VPD_LE{k}(:,3);
    Rnet_GS(sel)=NaN;
    
    LE_GS=ET;
    Bo_GS=H./LE;
    EF_GS=LE./(H+LE);
    LE_GS(sel)=NaN;
    Bo_GS(sel)=NaN;
    EF_GS(sel)=NaN;
    Bo_GS_AVG(k)=nanmedian(Bo_GS);
%     hist(Bo_GS);
%     title(num2str(k));
%     pause(1.5)
%     close all
    Bo_GS_MAX(k)=nanmedian(Bo_GS(VPD_GS>2));
    clear ind
    ind=find(~isnan(VPD_GS)&~isnan(LE_GS));
    
    x=VPD_GS(ind);
    y=LE_GS(ind);
    z=Bo_GS(ind);
    w=EF_GS(ind);
    
     for l =1:length(PRC_VPD)+1;
        if l == 1;
            ind=find(x<=PRC_VPD(l));
        elseif l<=length(PRC_VPD);
            ind=find(x>PRC_VPD(l-1) & x<=PRC_VPD(l));
        else
            ind=find(x>PRC_VPD(l-1));
        end
        if length(ind)<3 | isempty(ind)
            LE_BINS(k,l)=NaN;
            Bo_BINS(k,l)=NaN;
            EF_BINS(k,l)=NaN;
        else
            LE_BINS(k,l)=nanmedian(y(ind));
            Bo_BINS(k,l)=nanmedian(z(ind));
            EF_BINS(k,l)=nanmedian(w(ind));
        end
     end

end

par_gs(par_gs==0)=NaN;
par_gs(RMSE'./GS_BINS_UP(:,11)>0.1,:)=NaN;
ET_BINS_UP(ET_BINS_UP==0)=NaN;
GS_BINS_UP(GS_BINS_UP==0)=NaN;
