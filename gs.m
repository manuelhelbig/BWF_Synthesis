function [gs] = gs(H,LE,VPD,TA,ga,presstrue,lamdatrue,press)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% created by M. Helbig (helbigm@mcmaster.ca) on 2020-02-04
% calculates surface conductance (m s-1) by inverting the Penman-Monteith equation
% used for Helbig et al. (????) Increasing contribution of peatlands to boreal evapotranspiration in a warming climate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
% H: sensible heat flux (W m-2)
% LE: latent heat flux (W m-2)
% VPD: vapour pressure deficit (kPa)
% TA: air temperature (C)
% ga: aerodynamic conductance (m s-1)
% presstrue: 1 if measured barometric pressure should be used
% lamdatrue: 1 if latent heat of vaporization should be calculated based on TA
% press: barometric pressure (kPa)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% gas constant of dry air (J kg-1 K-1) (Foken 2008 p. 245)
Rd = 287.0586;        

% specific heat of air for constant pressure (J K-1 kg-1)
cp = 1004.834;  

% press: barometric pressure (kPa)
if presstrue==0;
    press=101.3;
end

% lambda: latent heat of vaporization (J kg-1)
if lamdatrue==1;
    lamda = (2.501 - 0.00237.*TA).*10^6;
else
    lamda=2.501.*10^6;
end

% air density (kg m-3)
rho = (press.*1000)./(Rd.*(TA+273.15));

% psychrometric constant in kPa K-1 
gamma = cp.*press./(0.622.*lamda);

% saturated vapour pressure (kPa) (Bolton 1980)
Esat = 611.2.*exp((17.62.*TA)./(243.12+TA))./1000;

% s -- slope of the saturation vapour pressure curve against temperature
% (kPa C-1)
s=(4098.*(0.6108.*exp((17.27.*TA)./(TA+237.3))))./(TA+237.4).^2;

% available energy (either H+LE or Rn-G)
% (W m-2)
Ra=H+LE;

% Ryu et al 2008 (equation 1)
eps = s./gamma;

% derivation of surface conductance (inversion of Penman-Monteith equation)
% e.g., Humphreys et al., 2007 (equation 2)
gs_inv = (ga).^(-1).*(eps.*(Ra./LE-1)-1)+((rho.*cp.*VPD)./(gamma.*LE));
gs = 1./gs_inv;
