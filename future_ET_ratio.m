%% estimating peatland and forest ET ratio %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% created by M. Helbig (helbigm@mcmaster.ca) on 2020-02-04
% used for Helbig et al. (????) Increasing contribution of peatlands to boreal evapotranspiration in a warming climate
function [ET_PTL_FUT45, ET_PTL_FUT85, ET_PTL_CUR, ET_FOR_FUT45, ET_FOR_FUT85, ET_FOR_CUR]...
    = future_ET_ratio(par_PTL,par_FOR,ga_ptl,ga_for,par_PAR_PTL,par_PAR_FOR,par_LAI,VPD_MAX_FUT45, VPD_MAX_FUT85, VPD_MAX_REC, dAV_45, dAV_85, LAI, RMSE_FOR, AV_JUL_DAY, FOR_AV)
%%% input variables %%%
% VPD_MAX_FUT45 = ESM projections (RCP4.5) of maximum VPD for 8 ESMs (structure with 8 layers)
% VPD_MAX_FUT85 = ESM projections (RCP8.5) of maximum VPD for 8 ESMs (structure with 8 layers)
% dAV_45 = change in mean daily available energy (H+LE) between 2006-2015 and 2091-2100 from 8 ESMs (RCP4.5)
% dAV_85 = change in mean daily available energy (H+LE) between 2006-2015 and 2091-2100 from 8 ESMs (RCP8.5)
% par_PTL = gs model parameters for peatlands
% par_FOR = gs model parameters for forests
% ga_ptl = aerodynamic conductance estimates for peatlands
% ga_for = aerodynamic conductance estimates for forests
% par_LAI = parameters for relationship between LAI and forest gs-VPD parameters
% par_PAR_FOR = relationship between g1 and g0 parameter (gs-VPD model) for forests
% par_PAR_PTL = relationship between g1 and g0 parameter (gs-VPD model) for peatlands
% RMSE_FOR = root mean squared error of gs model paramter vs LAI function
% AV_JUL_DAY = mean daily available energy for July (2004-2013 from FLUXCOM)
% FOR_AV = regression parameters of mean daily vs. mean afternoon available energy for forests
% PTL_AV = regression parameters of mean daily vs. mean afternoon available energy for forests

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% input variables %%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% boundary conditions
Ta = 20;            % air temperature, degC
cp = 1004.834;      % Specific heat of air, J kg-1 K-1        
Rstar = 8.3144;     % universal gas constant, J mol-1 K-1
Mair = 29;          % molecular mass of air, g mol-1  
pa = 1013.*100;     % atmospheric pressure (Pa)

lmbda = 3149000 - 2370 * (Ta+273.15);          % latent heat of vaporization
psych = (1013.*100 * cp) / (0.622 * lmbda);    % pyschrometric constant (Pa K-1)

air_density = 101.3 * Mair / (Rstar * (Ta+273.15)); % air density

% Slope of the saturation vapor pressure-temperature curve
% Thermodynamic relationship from HESS    (Pa K-1)
esat=100.*exp(54.8781919 - 6790.4985./ (Ta+273.15) - 5.02808.* log(Ta+273.15));

slope = esat * lmbda * 18./ (Rstar * (Ta+273.15) * (Ta+273.15) * 1000);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gs response to VPD
gsFxn=@(params,VPD) (params(2)+(1+params(1)./(sqrt(VPD))));

% number of samples for uncertainty estimation
sample_n=120;

% use standard deviation of peatland gs model parameters for uncertainty estimation
STD_PTL=std(par_PTL(:,1));

r=0;
% use 8 VPD projections to account for intermodel (ESM) differences
for v=1:8
    % current VPD (2006-2015) [kPa]
    dVPD_REC=VPD_MAX_REC;
    % future VPD (2091-2100)[kPa]
    dVPD_FUT=VPD_MAX_FUT45(:,:,v);
    dVPD_FUT_85=VPD_MAX_FUT85(:,:,v);
    
    % sample from parameter distribution to derive uncertainties
    for l=1:sample_n;
            
                % index for parameter selection
                r=r+1;
                r
                sel1 = round(rand(1).*length(par_FOR(:,1)));
                sel2 = round(rand(1).*length(par_PTL(:,1)));
                if sel1==0;
                   sel1=sel1+1;
                end
                if sel2==0;
                   sel2=sel2+1;
                end   
                % select aerodynamic conductance (15h to 18h during growing season)           
                ga_FOR=ga_for(sel1);
                ga_PTL=ga_ptl(sel2);
                
                % peatland: gs-VPD response
                % randomly select from parameter space for g1
                par_ptl=normrnd(nanmedian(par_PTL(:,1)),STD_PTL);
                % apply relationship between g1 and g0
                par_ptl(2)=par_PAR_PTL(2)+par_PAR_PTL(1).*par_ptl(1);
                
                % choose grid cells with available VPD data
                [index1 index2] = find(~isnan(dVPD_FUT) & ~isnan(dVPD_REC) & sum(isnan(VPD_MAX_FUT45),3)==0);
                for w=1:length(index1);
                    k=index1(w);
                    q=index2(w);
                    
                    % scale gs model parameters using LAI (for forests)
                    LAI_sel=LAI(k,q);
                    
                    % derive gs-VPD parameter estimates from LAI relationship
                    % add random error estimate to LAI relationship
                    std_par_g1(1) = normrnd(0,RMSE_FOR);
                    par_for(1)=(par_LAI(2)+LAI_sel.*par_LAI(1))+std_par_g1;
                    par_for(2)=par_PAR_FOR(2)+par_PAR_FOR(1).*par_for(1);
                    
                    % peatland gs for VPD (current climate)
                    GS_PTL_REC=gsFxn(par_ptl,dVPD_REC(k,q));
                    
                    % forest gs for VPD (current climate)
                    GS_FOR_REC=gsFxn(par_for,dVPD_REC(k,q));
                    
                    % mean available energy (July) for current climate
                    AV_EN=(AV_JUL_DAY(k,q)).*FOR_AV(1)+FOR_AV(2);
                    
                    % calculate grid cell foresy ET using Penman-Monteith (2006-2015)
                    ET_FOR_CUR(k,q,r)=(slope.*AV_EN+air_density.*cp.*(dVPD_REC(k,q).*1000).*ga_FOR)...
                        ./(slope+psych.*(1+ga_FOR./(GS_FOR_REC./1000)))./lmbda.*60.*60;
                    
                    AV_EN=(AV_JUL_DAY(k,q)).*PTL_AV(1)+PTL_AV(2);
                    ET_PTL_CUR(k,q,r)=(slope.*AV_EN+air_density.*cp.*(dVPD_REC(k,q).*1000).*ga_PTL)...
                        ./(slope+psych.*(1+ga_PTL./(GS_PTL_REC./1000)))./lmbda.*60.*60;
                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RCP4.5

                    % peatland and forest gs for VPD (future, RCP4.5)
                    GS_FOR_FUT=gsFxn(par_for,dVPD_FUT(k,q)); 
                    GS_PTL_FUT=gsFxn(par_ptl,dVPD_FUT(k,q)); 
                    
                    % estimate future available energy for forest
                    AV_EN=(dAV_45(k,q,v).*AV_JUL_DAY(k,q)).*FOR_AV(1)+FOR_AV(2);
                    % calculate grid cell forest ET (2091-2100) using Penman-Monteith
                    ET_FOR_FUT45(k,q,r)=real((slope.*AV_EN+air_density.*cp.*(dVPD_FUT(k,q).*1000).*ga_FOR)...
                        ./(slope+psych.*(1+ga_FOR./(GS_FOR_FUT./1000)))./lmbda.*60.*60);
                    
                    AV_EN=(dAV_45(k,q,v).*AV_JUL_DAY(k,q)).*PTL_AV(1)+PTL_AV(2);
                    ET_PTL_FUT45(k,q,r)=real((slope.*AV_EN+air_density.*cp.*(dVPD_FUT(k,q).*1000).*ga_PTL)...
                        ./(slope+psych.*(1+ga_PTL./(GS_PTL_FUT./1000)))./lmbda.*60.*60);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RCP8.5
                    % peatland and forest gs for VPD (future, RCP8.5)
                    GS_FOR_FUT=gsFxn(par_for,dVPD_FUT_85(k,q)); 
                    GS_PTL_FUT=gsFxn(par_ptl,dVPD_FUT_85(k,q)); 
                    
                    % calculate grid cell forest ET using Penman-Monteith
                    AV_EN=(dAV_85(k,q,v).*AV_JUL_DAY(k,q)).*FOR_AV(1)+FOR_AV(2);
                    ET_FOR_FUT85(k,q,r)=real((slope.*AV_EN+air_density.*cp.*(dVPD_FUT_85(k,q).*1000).*ga_FOR)...
                        ./(slope+psych.*(1+ga_FOR./(GS_FOR_FUT./1000)))./lmbda.*60.*60);
                    
                    AV_EN=(dAV_85(k,q,v).*AV_JUL_DAY(k,q)).*PTL_AV(1)+PTL_AV(2);
                    ET_PTL_FUT85(k,q,r)=real((slope.*AV_EN+air_density.*cp.*(dVPD_FUT_85(k,q).*1000).*ga_PTL)...
                        ./(slope+psych.*(1+ga_PTL./(GS_PTL_FUT./1000)))./lmbda.*60.*60);
                    clear AV_EN 
                end
                
                [index1 index2] = find(isnan(dVPD_FUT) & isnan(dVPD_REC));
                for w=1:length(index1);
                    k=index1(w);
                    q=index2(w);
                    ET_FOR_CUR(k,q,r)=NaN;
                    ET_PTL_CUR(k,q,r)=NaN;
                    ET_FOR_FUT45(k,q,r)=NaN;
                    ET_PTL_FUT45(k,q,r)=NaN;
                end
                 [index1 index2] = find(isnan(dVPD_FUT_85) & isnan(dVPD_REC));
                 for w=1:length(index1);
                    k=index1(w);
                    q=index2(w);
                    ET_FOR_FUT85(k,q,r)=NaN;
                    ET_PTL_FUT85(k,q,r)=NaN;
                end
                clear sel1 sel2
                if r==sample_n1
                    r=0;
                end
    end
end
