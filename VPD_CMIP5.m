%% load CMIP5 data %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% created by M. Helbig (helbigm@mcmaster.ca) on 2020-02-04
% calculate vapour pressure deficit from CMIP5 output
% used for Helbig et al. (????) Increasing contribution of peatlands to boreal evapotranspiration in a warming climate
function [LAT_CMIP5, LON_CMIP5, VPD_CMIP5] = VPD_CMIP5(path_HUSS, path_TAMAX, path_PSL)
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input:
% path_HUSS: path to specific humidity data
% path_TAMAX: path to maximum air temperature data
% path_PSL: path to near-surface atmospheric pressure data
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get specific humidity
finfo = dir([path_HUSS '*.nc']);
clear lat_CMIP5 lon_CMIP5

for k = 1:length(finfo);
    % read file in folder
    FILENAME = [path_HUSS finfo(k).name];
    finfo(k).name
    % in first iteration: read lat/long info
    if k == 1;
        latData  = ncread(FILENAME,'lat');
        lonData  = ncread(FILENAME,'lon');
        lonData(lonData>180)=lonData(lonData>180)-360;
        for n = 1:length(lonData);
            LAT_CMIP5(n,:) = latData;
        end

        for n = 1:length(latData);
            LON_CMIP5(:,n) = lonData;
        end
    end
    % read near-surface specific humidity
    HUSS  = ncread(FILENAME,'huss');  
    
    % get data and stack it (create subsets if needed)
    if k == 1;       
        HUSS_tmp = HUSS;
    else
        HUSS_tmp = cat(3,HUSS_tmp,HUSS);
    end
    clear HUSS
end

HUSS_CMIP5 = HUSS_tmp;
clear HUSS_tmp

% get near-surface maximum air temperature
finfo = dir([path_TAMAX '*.nc']);

for k = 1:length(finfo);
    % read file in folder
    FILENAME = [path_TAMAX finfo(k).name];
    finfo(k).name
    
    % read maximum air temperature
    TASMAX  = ncread(FILENAME,'tasmax');  
   
    % get data and stack it
    if k == 1;       
        TASMAX_tmp = TASMAX;
    else
        TASMAX_tmp = cat(3,TASMAX_tmp,TASMAX);
    end
    clear TASMAX
end

TASMAX_CMIP5 = TASMAX_tmp-273.15;
clear TASMAX_tmp

% get near-surface air pressure
finfo = dir([path_PSL '*.nc']);

for k = 1:length(finfo);
    % read file in folder
    FILENAME = [path_PSL finfo(k).name];
    finfo(k).name
    
    % read near-surface pressure
    PSL  = ncread(FILENAME,'psl');  
   
    % get data and stack it
    if k == 1;       
        PSL_tmp = PSL;
    else
        PSL_tmp = cat(3,PSL_tmp,PSL);
    end
    clear PSL
end

PSL_CMIP5 = PSL_tmp;
clear PSL_tmp

% derive VPD from specific humidity, pressure, air temperature
for k = 1:length(HUSS_CMIP5(1,1,:));
    for n = 1:100;
        if n == 1;
            e(:,:,n) = (HUSS_CMIP5(:,:,k).*PSL_CMIP5(:,:,k))/0.622;
        else
            e(:,:,n) = (HUSS_CMIP5(:,:,k)/0.622).*(PSL_CMIP5(:,:,k)-(0.622.*e(:,:,n-1)));
        end
    end
    e_vap(:,:,k) = e(:,:,end);
end

for k = 1:length(TASMAX_CMIP5(1,1,:)); 
    es_max(:,:,k) = 611.*exp((17.27.*TASMAX_CMIP5(:,:,k))./(237.3+TASMAX_CMIP5(:,:,k)));
    VPD_CMIP5(:,:,k) = ((es_max(:,:,k)./1000)-(e_vap(:,:,k)./1000));
end
clear e_vap es
