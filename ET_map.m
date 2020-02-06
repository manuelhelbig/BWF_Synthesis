%% map ratio of peatland and forest ET %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% created by M. Helbig (helbigm@mcmaster.ca) on 2020-02-04
% used for Helbig et al. (????) Increasing contribution of peatlands to boreal evapotranspiration in a warming climate
function ET_map(path, lat,lon, ET_PTL, ET_FOR, ET_PTL_FUT, ET_FOR_FUT)
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input:
% path: path to shapefile of boreal biome ['CIRCUM_BOREAL.shp']
% lat: latitude grid
% lon: longitude grid
% ET_PTL: gridded peatland evapotranspiration estimate
% ET_FOR: gridded forest evapotranspiration estimate
% ET_PTL_FUT: gridded future peatland evapotranspiration estimate
% ET_FOR_FUT: gridded future forest evapotranspiration estimate
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot ratio of peatland to forest ET
figure,
subplot(1,1,1);
ax = worldmap([43 90],[-180 180]);
land = shaperead([path 'CIRCUM_BOREAL.shp'], 'UseGeoCoords', true);
geoshow(ax, land, 'FaceColor', [176/255 196/255 222/255],'EdgeColor',[25/255 25/255 112/255]);
geoshow('landareas.shp','DefaultFaceColor', 'none', ...
        'DefaultEdgeColor', 'black','LineWidth',1.5);
hold on
l2 = geoshow(lat,lon,ET_PTL./ET_FOR,'DisplayType','texturemap');
l2.AlphaDataMapping = 'none'; % interpet alpha values as transparency values
l2.FaceAlpha = 'texturemap'; % Indicate that the transparency can be different each pixel
alpha(l2,double(~isnan(ET_PTL./ET_FOR))) % Change transparency to matrix where if data==NaN --> transparency = 1, else 0.
h_bar = colorbar('location','SouthOutside','FontSize',14);
set(get(h_bar,'ylabel'),'String', {'ET_{PTL}/ET_{FOR}'},'FontSize',14);
caxis([0.75 1.25]) 
colormap(flipud(brewermap([10],'RdGy')))
% plot projected change in ratio
subplot(1,2,2);
ax = worldmap([43 90],[-180 180]);
geoshow(ax, land, 'FaceColor', [176/255 196/255 222/255],'EdgeColor',[25/255 25/255 112/255]);
geoshow('landareas.shp','DefaultFaceColor', 'none', ...
        'DefaultEdgeColor', 'black','LineWidth',1.5);
hold on
l2 = geoshow(lat,lon,ET_PTL_FUT./ET_FOR_FUT-ET_PTL./ET_FOR,'DisplayType','texturemap');
l2.AlphaDataMapping = 'none'; % interpet alpha values as transparency values
l2.FaceAlpha = 'texturemap'; % Indicate that the transparency can be different each pixel
alpha(l2,double(~isnan(ET_PTL_FUT./ET_FOR_FUT))) % Change transparency to matrix where if data==NaN --> transparency = 1, else 0.
h_bar = colorbar('location','SouthOutside','FontSize',14);
set(get(h_bar,'ylabel'),'String', {'\Delta ET_{PTL}/ET_{FOR}'},'FontSize',14);
caxis([0 0.2]) 
colormap(brewermap([10],'YlOrRd'))
