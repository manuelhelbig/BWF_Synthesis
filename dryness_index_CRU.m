%% plotting dryness index from CRU data %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% created by M. Helbig (helbigm@mcmaster.ca) on 2020-02-04
% uses Climate Research Unit data to calculate dryness index
% used for Helbig et al. (????) Increasing contribution of peatlands to boreal evapotranspiration in a warming climate
function [LAT_CRU, LON_CRU, DI_CRU] = dryness_index_CRU(path_PRE, path_PET)
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input:
% path_PRE: path to folder containing files with precipitation data
% path_PET: path to folder containing files with potential evapotranspiration data
% CRU TS data can be downloaded from the Climate Research Unit, University of East Anglia website: https://crudata.uea.ac.uk/cru/data/hrg/cru_ts_4.03/
%%%%%%%%%%%%%%%%%%%%%%%%%%%
finfo = dir([path_PRE '*.nc']);
PRE=[];
for k = 1:length(finfo);
    if k == 1;
        lat = ncread([path_PRE,finfo(k).name],'lat');
        lon = ncread([path_PRE,finfo(k).name],'lon');
    end
    pre=ncread([path_PRE,finfo(k).name],'pre');
    if k > 1;
        PRE(:,:,length(PRE(1,1,:))+1:length(PRE(1,1,:))+length(pre(1,1,:))) = pre;
    else
        PRE = pre;
    end
    clear pre
end

finfo = dir([path_PET '*.nc']);
PET=[];
for k = 1:length(finfo);
    pet=ncread([path_PET,finfo(k).name],'pet');
    if k > 1;
        PET(:,:,length(PET(1,1,:))+1:length(PET(1,1,:))+length(pet(1,1,:))) = pet;
    else
        PET = pet;
    end
    clear pet
end

MM=[];
YY=[];
start_yr=str2num(finfo(1).name(12:15));
for k=1:length(PRE(1,1,:))/12;
    MM=horzcat(MM,1:12);
    YY=horzcat(YY,ones(1,12).*(start_yr+k));
end

for k = 1:length(PRE(1,1,:))/12;
    PET_ANN(:,:,k)=nanmean(PET(:,:,YY==start_yr+k),3).*365;
    PRE_ANN(:,:,k)=nansum(PRE(:,:,YY==start_yr+k),3);
    PRE_GS(:,:,k)=nansum(PRE(:,:,YY==start_yr+k & (MM>=5 & MM<=9)),3);
    PET_GS(:,:,k)=nanmean(PET(:,:,YY==start_yr+k & (MM>=5 & MM<=9)),3).*153;
end
YY=start_yr:start_yr+length(PRE(1,1,:))/12-1;

lat = double(lat);
lon = double(lon);
for k = 1:720;
LAT_CRU(k,:) = lat;
end

for k = 1:360;
LON_CRU(:,k) = lon;
end

avg_start=1981;
avg_end=2010;
DI_CRU=nanmean(PET_ANN(:,:,YY>=avg_start & YY<=avg_end),3)./...
    nanmean(PRE_ANN(:,:,YY>=avg_start & YY<=avg_end),3);
