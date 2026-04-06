% The purpose of this script is to create correlation maps of CPUE and
% percent catch with oceanographic and climate index data.
% It borrows heavily from existing scripts.

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
Vessels = ncread('../FisheryData/TotalVessels.nc', 'Total Vessels');
Effort = ncread('../FisheryData/TotalEffort.nc', 'Total Effort');
Bigeye = ncread('../FisheryData/TotalBigeyeCaught.nc', 'Total Bigeye Caught');
Yellowfin = ncread('../FisheryData/TotalYellowfinCaught.nc', 'Total Yellowfin Caught');
Mahi = ncread('../FisheryData/TotalMahiCaught.nc', 'Total Mahi Caught');
Pomfret = ncread('../FisheryData/TotalPomfretCaught.nc', 'Total Pomfret Caught');
Sword = ncread('../FisheryData/TotalSwordfishCaught.nc', 'Total Swordfish Caught');
PMUS = ncread('../FisheryData/TotalPMUS.nc', 'Total PMUS Caught');
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
% Vessels = ncread('..\FisheryData\TotalVessels.nc', 'Total Vessels');
% Effort = ncread('..\FisheryData\TotalEffort.nc', 'Total Effort');
% Bigeye = ncread('..\FisheryData\TotalBigeyeCaught.nc', 'Total Bigeye Caught');
% Yellowfin = ncread('..\FisheryData\TotalYellowfinCaught.nc', 'Total Yellowfin Caught');
% Mahi = ncread('..\FisheryData\TotalMahiCaught.nc', 'Total Mahi Caught');
% Pomfret = ncread('..\FisheryData\TotalPomfretCaught.nc', 'Total Pomfret Caught');
% Sword = ncread('..\FisheryData\TotalSwordfishCaught.nc', 'Total Swordfish Caught');
% PMUS = ncread('..\FisheryData\TotalPMUS.nc', 'Total PMUS Caught');
% ONI = readtable('..\ClimateIndices\ONI_withPhases.csv');
% PDO = readtable('..\ClimateIndices\PDO.csv');
% NPGO = readtable('..\ClimateIndices\NPGO.csv');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Wrangle data
% Permute so arrays are Lat x Lon (x Depth) x time
Catchability = permute(Catchability, [2 1 3]);
Depth8 = permute(Depth8, [2 1 3]);
Depth14 = permute(Depth14, [2 1 3]);
Thickness814 = permute(Thickness814, [2 1 3]);
O2_2mlpl = permute(O2_2mlpl, [2 1 3]);
Vessels = permute(Vessels, [2 1 3]);
Effort = permute(Effort, [2 1 3]);
Bigeye = permute(Bigeye, [2 1 3]);
Yellowfin = permute(Yellowfin, [2 1 3]);
Mahi = permute(Mahi, [2 1 3]);
Pomfret = permute(Pomfret, [2 1 3]);
Sword = permute(Sword, [2 1 3]);
PMUS = permute(PMUS, [2 1 3]);

% Limit climate indices to our period of interest: 1995–2024
ONI_9524 = ONI(ONI.YR >= 1995 & ONI.YR <= 2024,:);
PDO_9524 = PDO(PDO.Year >= 1995 & PDO.Year <= 2024,:);
NPGO_9524 = NPGO(NPGO.YEAR >= 1995 & NPGO.YEAR <= 2024,:);
% Clean up
clear ONI PDO NPGO

% Create lat and lon grids (for plotting)
% Horizontally stack column with latitudes
lat_grid = repmat(Lat,1,size(Lon,1));
% Transpose column with longitudes to a row and vertically stack
lon_grid = repmat(Lon',size(Lat,1),1);

% Sum vessels in each grid cell over time
% To use as a filter for plotting
Vessels_total = sum(Vessels, 3, "omitnan");

% Calculate CPUEs and percent of PMUS
BigeyeCPUE = Bigeye ./ Effort * 1000;
PctBigeye = Bigeye ./ PMUS * 100;
YellowfinCPUE = Yellowfin ./ Effort * 1000;
PctYellowfin = Yellowfin ./ PMUS * 100;
MahiCPUE = Mahi ./ Effort * 1000;
PctMahi = Mahi ./ PMUS * 100;
PomfretCPUE = Pomfret ./ Effort * 1000;
PctPomfret = Pomfret ./ PMUS * 100;
SwordCPUE = Sword ./ Effort * 1000;
PctSword = Sword ./ PMUS * 100;
PMUSCPUE = PMUS ./ Effort * 1000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Create correlation maps
% This time around, putting the correlation step into the plotting step to
% avoid making a bunch of empty matrices to fill.  Might be slower, but
% also might be less prone to error and faster to write the code for.

figure
Corr_Map(lat_grid, lon_grid, PDO_9524.PDO, BigeyeCPUE, Vessels_total, 0.05, ...
    -1, 1, 'PDO and Bigeye CPUE')
% saveas(gcf, 'PDO_BigeyeCPUE_Pearson.pdf')

figure
Corr_Map(lat_grid, lon_grid, PDO_9524.PDO, PomfretCPUE, Vessels_total, 0.05, ...
    -1, 1, 'PDO and Pomfret CPUE')
% saveas(gcf, 'PDO_PomfretCPUE_Pearson.pdf')

figure
Corr_Map(lat_grid, lon_grid, PDO_9524.PDO, MahiCPUE, Vessels_total, 0.05, ...
    -1, 1, 'PDO and Mahi CPUE')
% saveas(gcf, 'PDO_MahiCPUE_Pearson.pdf')

figure
Corr_Map(lat_grid, lon_grid, PDO_9524.PDO, YellowfinCPUE, Vessels_total, 0.05, ...
    -1, 1, 'PDO and Yellowfin CPUE')
saveas(gcf, 'PDO_YellowfinCPUE_Pearson.pdf')

figure
Corr_Map(lat_grid, lon_grid, NPGO_9524.NPGO, BigeyeCPUE, Vessels_total, 0.05, ...
    -1, 1, 'NPGO and Bigeye CPUE')
% saveas(gcf, 'NPGO_BigeyeCPUE_Pearson.pdf')

figure
Corr_Map(lat_grid, lon_grid, NPGO_9524.NPGO, YellowfinCPUE, Vessels_total, 0.05, ...
    -1, 1, 'NPGO and Yellowfin CPUE')
% saveas(gcf, 'NPGO_YellowfinCPUE_Pearson.pdf')

figure
Corr_Map(lat_grid, lon_grid, NPGO_9524.NPGO, MahiCPUE, Vessels_total, 0.05, ...
    -1, 1, 'NPGO and Mahi CPUE')
% saveas(gcf, 'NPGO_MahiCPUE_Pearson.pdf')

figure
Corr_Map(lat_grid, lon_grid, NPGO_9524.NPGO, PomfretCPUE, Vessels_total, 0.05, ...
    -1, 1, 'NPGO and Pomfret CPUE')
saveas(gcf, 'NPGO_PomfretCPUE_Pearson.pdf')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Functions
% Function to calculate correlations and plot significant ones
function Corr_Map(lat_grid_used, lon_grid_used, env_var, catch_var, vessels, ...
    sig_level, min_val, max_val, map_title)

% Empty matrices to fill
r_vals(1:size(lat_grid_used,1), 1:size(lon_grid_used,2)) = NaN;
p_vals(1:size(lat_grid_used,1), 1:size(lon_grid_used,2)) = NaN;

% Pearson correlation, loop through grid cells
for r = 1:1:size(lat_grid_used,1)
    for c = 1:1:size(lon_grid_used,2)

        % Climate indices
        if size(env_var, 3) == 1
            [C_r, C_p] = corrcoef(env_var(:,1), catch_var(r,c,:), 'Rows', 'pairwise');
            r_vals(r, c) = C_r(1,2);
            p_vals(r, c) = C_p(1,2);
        end

        % Oceanographic variables
        if size(env_var, 3) == 360
            [C_r, C_p] = corrcoef(env_var(r,c,:), catch_var(r,c,:), 'Rows', 'pairwise');
            r_vals(r, c) = C_r(1,2);
            p_vals(r, c) = C_p(1,2);
        end
    end
end

% Only plot significant values
r_vals(p_vals >= sig_level) = NaN;

% Omit confidential cells
r_vals(vessels < 3) = NaN;

% Map
axesm('mercator','MapLatLimit',[10 40],'MapLonLimit',[180 230], ...
    'MLineLocation', 10, 'PLineLocation', 10, ... % draw every 10 degrees         
    'Grid', 'on', 'MeridianLabel','on','ParallelLabel','on', ...
    'MLabelParallel', 10, 'MLabelLocation', 10, 'PLabelLocation', 10); % label every 10 degrees, below map 
Ig = geoshow(lat_grid_used,lon_grid_used,r_vals,'DisplayType','texturemap');
set(Ig,'AlphaData',double(~isnan(r_vals)),'AlphaDataMapping','none','FaceAlpha','texturemap'); 
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