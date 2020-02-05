function [ga] = ga(USTAR,WS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% created by M. Helbig (helbigm@mcmaster.ca) on 2020-02-04
% calculates surface conductance (m s-1) by inverting the Penman-Monteith equation
% used for Helbig et al. (????) Increasing contribution of peatlands to boreal evapotranspiration in a warming climate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
% USTAR: friction velocity (m s-1)
% WS: wind speed (m s-1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% see for example Humphreys et al., 2007 (equation 3)
ga =((2./(0.4.*USTAR)).*(0.89).^(2/3)+(WS./(USTAR.^2))).^-1;
