% The purpose of this script is to look at composite oceanographic states
% for each phase of each mode of climate variability.
% We'll also look at spatiatemploral correlations because much of the prep
% work is the same.
% Update: removed some of the earlier code (composite states) and
% streamlined what remains.  I was going to post a new script but then I
% decided that this is what version control is for, right?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Load data
% Matlab needs a "here" package.  Until then, use the path that fits your
% OS.
% Mac:
Lat = ncread('../OceanData/BigeyeCatchability.nc', 'Latitude');
Lon = ncread('../OceanData/BigeyeCatchability.nc', 'Longitude');
Month = ncread('../OceanData/BigeyeCatchability.nc', 'Month');
Catchability = ncread('../OceanData/BigeyeCatchability.nc', 'Catchability');
Depth8 = ncread('../OceanData/Depth8deg.nc', 'Depth8deg');
Depth14 = ncread('../OceanData/Depth14deg.nc', 'Depth14deg');
Thickness814 = ncread('../OceanData/Thickness814deg.nc', 'Thickness814deg');
O2_2mlpl = ncread('../OceanData/O2_2mlpl_depth_1deg_noInterp.nc', 'GLORYS_O2_2MLPL_REGRID');
ONI = readtable('../ClimateIndices/ONI_withPhases.csv');
PDO = readtable('../ClimateIndices/PDO.csv');
NPGO = readtable('../ClimateIndices/NPGO.csv');
% PC:
% Lat = ncread('..\OceanData\BigeyeCatchability.nc', 'Latitude');
% Lon = ncread('..\OceanData\BigeyeCatchability.nc', 'Longitude');
% Month = ncread('..\OceanData\BigeyeCatchability.nc', 'Month');
% Catchability = ncread('..\OceanData\BigeyeCatchability.nc', 'Catchability');
% Depth8 = ncread('..\OceanData\Depth8deg.nc', 'Depth8deg');
% Depth14 = ncread('..\OceanData\Depth14deg.nc', 'Depth14deg');
% Thickness814 = ncread('..\OceanData\Thickness814deg.nc', 'Thickness814deg');
% O2_2mlpl = ncread('..\OceanData\O2_2mlpl_depth_1deg_noInterp.nc', 'GLORYS_O2_2MLPL_REGRID');
% ONI = readtable('..\ClimateIndices\ONI_withPhases.csv');
% PDO = readtable('..\ClimateIndices\PDO.csv');
% NPGO = readtable('..\ClimateIndices\NPGO.csv');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Wrangle data
% Limit climate indices to our period of interest: 1995–2024
ONI_9524 = ONI(ONI.YR >= 1995 & ONI.YR <= 2024,:);
PDO_9524 = PDO(PDO.Year >= 1995 & PDO.Year <= 2024,:);
NPGO_9524 = NPGO(NPGO.YEAR >= 1995 & NPGO.YEAR <= 2024,:);
% Clean up
clear ONI PDO NPGO

% Permute so arrays are Lat x Lon x time
Catchability = permute(Catchability, [2 1 3]);
Depth8 = permute(Depth8, [2 1 3]);
Depth14 = permute(Depth14, [2 1 3]);
Thickness814 = permute(Thickness814, [2 1 3]);
O2_2mlpl = permute(O2_2mlpl, [2 1 3]);

% Create lat and lon grids (for plotting)
% Horizontally stack column with latitudes
lat_grid = repmat(Lat,1,size(Lon,1));
% Transpose column with longitudes to a row and vertically stack
lon_grid = repmat(Lon',size(Lat,1),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Create correlation maps
% This time around, putting the correlation step into the plotting step to
% avoid making a bunch of empty matrices to fill.  Might be slower, but
% also might be less prone to error and faster to write the code for.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plots
% There are functions for this at the end of the script
figure
Corr_Ctr_Map(lat_grid, lon_grid, PDO_9524.PDO, O2_2mlpl, O2_2mlpl, 100, 100, ...
    700, 0.05, -1, 1, 'PDO & O_2')
saveas(gcf, 'PDO_O2_Pearson_Contoured.pdf')

figure
Corr_Ctr_Map(lat_grid, lon_grid, PDO_9524.PDO, Depth8, Depth8, 100, 50, ...
    700, 0.05, -1, 1, 'PDO & Depth 8-deg')
saveas(gcf, 'PDO_D8_Pearson_Contoured.pdf')

figure
Corr_Ctr_Map(lat_grid, lon_grid, PDO_9524.PDO, Depth14, Depth14, 100, 50, ...
    700, 0.05, -1, 1, 'PDO & Depth 14-deg')
saveas(gcf, 'PDO_D14_Pearson_Contoured.pdf')

figure
Corr_Ctr_Map(lat_grid, lon_grid, NPGO_9524.NPGO, O2_2mlpl, O2_2mlpl, 100, 100, ...
    700, 0.05, -1, 1, 'NPGO & O_2')
saveas(gcf, 'NPGO_O2_Pearson_Contoured.pdf')

figure
Corr_Ctr_Map(lat_grid, lon_grid, NPGO_9524.NPGO, Depth8, Depth8, 100, 50, ...
    700, 0.05, -1, 1, 'NPGO & Depth 8-deg')
saveas(gcf, 'NPGO_D8_Pearson_Contoured.pdf')

figure
Corr_Ctr_Map(lat_grid, lon_grid, NPGO_9524.NPGO, Depth14, Depth14, 100, 50, ...
    700, 0.05, -1, 1, 'NPGO & Depth 14-deg')
saveas(gcf, 'NPGO_D14_Pearson_Contoured.pdf')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Functions
% Function to calculate correlations and plot significant ones, plus
% overlay contours
function Corr_Ctr_Map(lat_grid_used, lon_grid_used, clim_var, ocn_var, ...
    ctr_var, ctr_min, ctr_step, ctr_max, sig_level, min_val, max_val, map_title)

% Empty matrices to fill
r_vals(1:size(lat_grid_used,1), 1:size(lon_grid_used,2)) = NaN;
p_vals(1:size(lat_grid_used,1), 1:size(lon_grid_used,2)) = NaN;

% Pearson correlation, loop through grid cells
for r = 1:1:size(lat_grid_used,1)
    for c = 1:1:size(lon_grid_used,2)

        [C_r, C_p] = corrcoef(clim_var(:,1), ocn_var(r,c,:), 'Rows', 'pairwise');
        r_vals(r, c) = C_r(1,2);
        p_vals(r, c) = C_p(1,2);

    end
end

% Climatology for contouring
ctr_climo = mean(ctr_var, 3, "omitnan");

% Only plot significant values
r_vals(p_vals >= sig_level) = NaN;

% Map
axesm('mercator','MapLatLimit',[10 40],'MapLonLimit',[180 230], ...
    'MLineLocation', 10, 'PLineLocation', 10, ... % draw every 10 degrees         
    'Grid', 'on', 'MeridianLabel','on','ParallelLabel','on', ...
    'MLabelParallel', 10, 'MLabelLocation', 10, 'PLabelLocation', 10); % label every 10 degrees, below map 
Ig = geoshow(lat_grid_used,lon_grid_used,r_vals,'DisplayType','texturemap');
set(Ig,'AlphaData',double(~isnan(r_vals)),'AlphaDataMapping','none','FaceAlpha','texturemap'); 
[C,h] = contourm(lat_grid_used,lon_grid_used,ctr_climo, ctr_min:ctr_step:ctr_max, 'k');
clabelm(C,h)
% plotm([20 20], [180 230], 'k');
% plotm([10 40], [210 210], 'k');
% plotm([26 26], [180 210], 'k');
geoshow('landareas.shp','FaceColor',[0.5 0.5 0.5]);
colormap(redblueTecplot)
clim([min_val max_val])
cb = colorbar();
ylabel(cb, 'Correlation Coefficient', 'Rotation', 270);
title(sprintf('%s'), map_title)
set(gcf,'renderer','Painters')
tightmap
end

