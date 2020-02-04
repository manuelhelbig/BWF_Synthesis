%% boundary line analysis %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% created by M. Helbig (helbigm@mcmaster.ca) on 2020-02-04
% used for Helbig et al. (????) Increasing contribution of peatlands to boreal evapotranspiration in a warming climate
function [par_gs, RMSE_GS, R2_GS, VPD_BINS_UP, GS_BINS_UP, ET_BINS_UP, Bo_BINS, EF_BINS] = boundary_line_analysis(DOY, HH, H, LE, VPD, RN, SWIN, TA, GS, ga, time_start, time_end, DOY_start, DOY_end)
%%% input variables %%%
% DOY = day of year of dataset
% HH = timestamp (hours)
% H = sensible heat flux (W m-2)
% LE = latent heat flux (W m-2)
% VPD = vapour pressure deficit (kPa)
% RN = net radiation (W m-2)
% SWIN = incoming shortwave radiation (W m-2)
% TA = air temperature (degC)
% GS = surface conductance (m s-1)
% ga = aerodynamic conductance (m s-1)
% define time window (hours) for analysis
    % time_start = 15.5;
    % time_end = 18;
% define time window (day of year) for analysis
    % DOY_start=121;
    % DOY_end=273;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% input variables %%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% heat of vaporization [J kg-1]
Lv         = 2.5e6;                
opts = statset('nlinfit');
opts.RobustWgtFun = 'bisquare';
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
% convert ET to mm hour-1
ET = LE.*30.*60./Lv.*2;

% remove data for outside the window of analysis
VPD(sel1)=NaN;
GS(sel)=NaN;
RN(sel)=NaN;
AV(sel)=NaN;
ET(sel2)=NaN;
SWIN(sel)=NaN;
TA(sel)=NaN;
ga(sel)=NaN;
% get median aerodynamic conductance per site
ga_SIT=nanmedian(ga);

% VPD bins (increments of 0.1 kPa)
PRC_VPD=0.1:0.1:3;
i=find(~isnan(GS));

VPD_UP = [];
GS_UP = [];
Rn_UP = [];
ET_UP = [];
for l =1:length(PRC_VPD)+1;
        % select data within VPD bins
        if l == 1;
            ind=find(VPD<=PRC_VPD(l));     
        elseif l<=length(PRC_VPD);
            ind=find(VPD>PRC_VPD(l-1) & VPD<=PRC_VPD(l));  
        else
            ind=find(VPD>PRC_VPD(l-1));
        end

        if length(ind)<3 | isempty(i)
            % do not use surface conductance if only few datapoints are
            % available
            GS_BINS(l)=NaN;
            GS_BINS_SD(l)=NaN;  
            y1=NaN;

            VPD_sel = VPD(ind);
            RN_sel = RN(ind);
            AV_sel = AV(ind);
            ET_sel = ET(ind);

            % detect outliers using quartile method
            IQR = prctile(ET_sel,75)-prctile(ET_sel,25);
            outl = find(ET_sel>prctile(ET_sel,75)+1.5*IQR | ET_sel<prctile(ET_sel,25)-1.5*IQR);
            ET_sel(outl)=NaN;

            VPD_BINS(l)=nanmean(VPD_sel);
            RN_BINS(l)=nanmean(RN_sel);
            ET_BINS(l)=nanmean(ET_sel);
            ET_BINS_SD(l)=nanstd(ET_sel);
            AV_BINS(l)=nanmean(AV_sel);
            VPD_BINS_SD(l)=nanstd(VPD_sel);
            
            % define upper boundary (missing gs)
            x1=NaN;
            y1=NaN;
            z1=NaN;
            w1=NaN;
            
            VPD_BINS_UP(l)=NaN;
            GS_BINS_UP(l)=NaN;
            ET_BINS_UP(l)=NaN;
            AV_BINS_UP(l)=NaN;
            RN_BINS_UP(l)=NaN;
            
        else
            VPD_sel = VPD(ind);
            RN_sel = RN(ind);
            GS_sel = GS(ind);
            AV_sel = AV(ind);
            ET_sel = ET(ind);

            % detect outliers using quartile method
            IQR = prctile(GS_sel,75)-prctile(GS_sel,25);
            outl = find(GS_sel>prctile(GS_sel,75)+1.5*IQR | GS_sel<prctile(GS_sel,25)-1.5*IQR);
            GS_sel(outl)=NaN;

            % detect outliers using quartile method
            IQR = prctile(ET_sel,75)-prctile(ET_sel,25);
            outl = find(ET_sel>prctile(ET_sel,75)+1.5*IQR | ET_sel<prctile(ET_sel,25)-1.5*IQR);
            ET_sel(outl)=NaN;

            % mean per VPD bin
            VPD_BINS(l)=nanmean(VPD_sel);
            GS_BINS(l)=nanmean(GS_sel);
            ET_BINS(l)=nanmean(ET_sel);
            RN_BINS(l)=nanmean(RN_sel);
            AV_BINS(l)=nanmean(AV_sel);

            % standard deviation per VPD bin
            VPD_BINS_SD(l)=nanstd(VPD_sel);
            GS_BINS_SD(l)=nanstd(GS_sel);
            ET_BINS_SD(l)=nanstd(ET_sel);

            % define upper boundary (> mean + 1 std)
            x1=VPD_sel(GS_sel>=GS_BINS(l)+1.*GS_BINS_SD(l));
            y1=GS_sel(GS_sel>=GS_BINS(l)+1.*GS_BINS_SD(l));
            z1=ET_sel(GS_sel>=GS_BINS(l)+1.*GS_BINS_SD(l));
            w1=RN_sel(GS_sel>=GS_BINS(l)+1.*GS_BINS_SD(l));
            v1=AV_sel(GS_sel>=GS_BINS(l)+1.*GS_BINS_SD(l));
            clear GS_sel VPD_sel RN_sel ET_sel AV_sel

            % median of upper boundary per VPD bin
            VPD_BINS_UP(l)=nanmedian(x1);
            GS_BINS_UP(l)=nanmedian(y1);
            ET_BINS_UP(l)=nanmedian(z1);
            AV_BINS_UP(l)=nanmedian(v1);
            RN_BINS_UP(l)=nanmedian(w1);
        end
        % concatenate upper boundary data point
        VPD_UP=vertcat(VPD_UP,x1);
        GS_UP=vertcat(GS_UP,y1);
        Rn_UP=vertcat(Rn_UP,w1);
        ET_UP=vertcat(ET_UP,z1);
        clear x1 y1 w1 z1

end 
% fit model for VPD-gs relationship
        % only fit if more than 10 datapoint are available
        if sum(~isnan(GS_UP))>10;
            % only fit model to VPD > 1kPa
            thr=1;
            GS_UP(VPD_UP<thr)=NaN;
            VPD_UP(VPD_UP<thr)=NaN;
            Rn_UP(VPD_UP<thr)=NaN;

            x = VPD_UP;
            % convert to mm s-1
            y = GS_UP.*1000;

            i1 = find(~isnan(x) & ~isnan(y));
            % fit VPD-gs function (see eqn. ? in Helbig et al., ????)
            gsFxn=@(params,VPD) (params(2)+(1+params(1)./(sqrt(VPD))));
            % initial parameter estimates
            par0_gs=[10 -5];
            % non-linear model fit
            par_gs = nlinfit(x(i1),y(i1),gsFxn,par0_gs,opts);
            clear x y
            % goodness of fit
            x=VPD_BINS_UP;
            y=GS_BINS_UP;
            RMSE_GS=sqrt(mean(y(~isnan(x) & ~isnan(y) & x>thr)-(gsFxn(par_gs,x(~isnan(x) &~isnan(y) & x>thr))./1000)).^2);
            R=corrcoef(y(~isnan(x) & ~isnan(y) & x>thr),(gsFxn(par_gs,x(~isnan(x) &~isnan(y) & x>thr))./1000));
            R2_GS=R(1,2)^2;

        else
            RMSE_GS=NaN;
            par_gs = [NaN NaN];
            R2_GS=NaN;
        end    
        clear GS_UP VPD_UP

% get Bowen ratio
% select data for analysis window
sel = find((DOY<DOY_start | DOY>=DOY_end) | LE<50 ...
        | (HH<time_start | HH>=time_end) | (H+LE)<100);

VPD_GS=VPD;
LE_GS=ET;
% Bowen ratio
Bo_GS=H./LE;
% evaporative fraction
EF_GS=LE./(H+LE);

LE_GS(sel)=NaN;
Bo_GS(sel)=NaN;
EF_GS(sel)=NaN;
VPD_GS(sel)=NaN;
clear sel

ind=find(~isnan(VPD_GS)&~isnan(LE_GS));
x=VPD_GS(ind);
y=LE_GS(ind);
z=Bo_GS(ind);
w=EF_GS(ind);
% derive Bowen ratio response to VPD
 for l =1:length(PRC_VPD)+1;
    if l == 1;
        ind=find(x<=PRC_VPD(l));
    elseif l<=length(PRC_VPD);
        ind=find(x>PRC_VPD(l-1) & x<=PRC_VPD(l));
    else
        ind=find(x>PRC_VPD(l-1));
    end
    if length(ind)<3 | isempty(ind)
        LE_BINS(l)=NaN;
        Bo_BINS(l)=NaN;
        EF_BINS(l)=NaN;
    else
        LE_BINS(l)=nanmedian(y(ind));
        Bo_BINS(l)=nanmedian(z(ind));
        EF_BINS(l)=nanmedian(w(ind));
    end
 end
