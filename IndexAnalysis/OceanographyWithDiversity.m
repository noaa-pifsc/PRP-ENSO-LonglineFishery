% The purpose of this script is to look at correlations between measures of 
% species richness/diversity and oceanography.  
% Much of the code comes from CompositeOceanography.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Load data
% Matlab needs a "here" package.  Until then, use the path that fits your
% OS.
% Mac:
Lat = ncread('../OceanData/BigeyeCatchability.nc', 'Latitude');
Lon = ncread('../OceanData/BigeyeCatchability.nc', 'Longitude');
Month = ncread('../OceanData/BigeyeCatchability.nc', 'Month');
Catchability = ncread('../OceanData/BigeyeCatchability.nc', 'Catchability');
GODAS_Depth = ncread('../OceanData/GODAS_1deg.nc', 'LEV1_31');
GODAS = ncread('../OceanData/GODAS_1deg.nc', 'GODAS_REGRID');
O2_2mlpl = ncread('../OceanData/O2_2mlpl_depth_1deg_noInterp.nc', 'GLORYS_O2_2MLPL_REGRID');
Richness = ncread('../FisheryData/SpeciesRichness.nc', 'Species Richness');
Shannon = ncread('../FisheryData/ShannonIndex.nc', 'Shannon Index');
Simpson = ncread('../FisheryData/SimpsonIndex.nc', 'Simpson Index');
Vessels = ncread('../FisheryData/TotalVessels.nc', 'Total Vessels');
Effort = ncread('../FisheryData/TotalEffort.nc', 'Total Effort');
% PC:
% Lat = ncread('..\OceanData\BigeyeCatchability.nc', 'Latitude');
% Lon = ncread('..\OceanData\BigeyeCatchability.nc', 'Longitude');
% Month = ncread('..\OceanData\BigeyeCatchability.nc', 'Month');
% Catchability = ncread('..\OceanData\BigeyeCatchability.nc', 'Catchability');
% GODAS = ncread('..\OceanData\GODAS_1deg.nc', 'GODAS_REGRID');
% GODAS_Depth = ncread('..\OceanData\GODAS_1deg.nc', 'LEV1_31');
% O2_2mlpl = ncread('..\OceanData\O2_2mlpl_depth_1deg_noInterp.nc', 'GLORYS_O2_2MLPL_REGRID');
% ONI = readtable('..\ClimateIndices\ONI_withPhases.csv');
% PDO = readtable('..\ClimateIndices\PDO.csv');
% NPGO = readtable('..\ClimateIndices\NPGO.csv');
% Richness = ncread('..\FisheryData\SpeciesRichness.nc', 'Species Richness');
% Shannon = ncread('..\FisheryData\ShannonIndex.nc', 'Shannon Index');
% Simpson = ncread('..\FisheryData\SimpsonIndex.nc', 'Simpson Index');
% Vessels = ncread('..\FisheryData\TotalVessels.nc', 'Total Vessels');
% Effort = ncread('..\FisheryData\TotalEffort.nc', 'Total Effort');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Wrangle data
% Permute so arrays are Lat x Lon (x Depth) x time
Catchability = permute(Catchability, [2 1 3]);
GODAS = permute(GODAS, [2 1 3 4]);
O2_2mlpl = permute(O2_2mlpl, [2 1 3]);
Richness = permute(Richness, [2 1 3]);
Shannon = permute(Shannon, [2 1 3]);
Simpson = permute(Simpson, [2 1 3]);
Vessels = permute(Vessels, [2 1 3]);
Effort = permute(Effort, [2 1 3]);

% Create lat and lon grids (for plotting)
% Horizontally stack column with latitudes
lat_grid = repmat(Lat,1,size(Lon,1));
% Transpose column with longitudes to a row and vertically stack
lon_grid = repmat(Lon',size(Lat,1),1);

% Sum vessels in each grid cell over time
% To use as a filter for plotting
Vessels_total = sum(Vessels, 3, "omitnan");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Climatologies 
Catch_climo = mean(Catchability, 3, "omitnan");
O2_climo = mean(O2_2mlpl, 3, "omitnan");

Richness_climo = mean(Richness, 3, "omitnan");
Shannon_climo = mean(Shannon, 3, "omitnan");
Simpson_climo = mean(Simpson, 3, "omitnan");

Effort_total = sum(Effort, 3, "omitnan");

RichnessPUE = Richness ./ Effort * 10^6;

% One approach to temperature:
% Max depth of 8-deg C, min depth of 14-deg waters, and thickness of 8-14-deg C layer
% Empty matrices to fill
GODAS_8degDepth_climo(1:length(Lat), 1:length(Lon)) = NaN;
GODAS_14degDepth_climo(1:length(Lat), 1:length(Lon)) = NaN;
GODAS_8to14degThickness_climo(1:length(Lat), 1:length(Lon)) = NaN;

for r = 1:1:length(Lat)
    for c = 1:1:length(Lon)
        % Each grid cell has 31 depths at each of 360 months
        % At each time step, we're going to:
        % Find the max depth of 8-deg waters,
        % min depth of 14-deg waters, and
        % Determine the thickness of the 8-14-deg layer
        % Then, we'll average over time and fill the matrix
        Depth8(1:360) = NaN;
        Depth14(1:360) = NaN;
        Thick814(1:360) = NaN;
        for m = 1:1:360
            Z8_idx = find(GODAS(r,c,:,m) >= 8);
            Z14_idx = find(GODAS(r,c,:,m) <= 14);
            
            Depth8(m) = GODAS_Depth(max(Z8_idx));
            Depth14(m) = GODAS_Depth(min(Z14_idx));
            Thick814(m) = GODAS_Depth(max(Z8_idx)) - GODAS_Depth(min(Z14_idx));
            
        end
        clear m *idx

        GODAS_8degDepth_climo(r,c) = mean(Depth8, "omitnan");
        GODAS_14degDepth_climo(r,c) = mean(Depth14, "omitnan");
        GODAS_8to14degThickness_climo(r,c) = mean(Thick814, "omitnan");

    end
end
clear r c Depth8 Depth14 Thick814

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Correlations with diversity metrics
% Sticking to Pearson only following group decision
% Empty arrays to fill 
Catch_Richness_r(1:length(Lat), 1:length(Lon)) = NaN;
Catch_RichnessPUE_r(1:length(Lat), 1:length(Lon)) = NaN;
Catch_Shannon_r(1:length(Lat), 1:length(Lon)) = NaN;
Catch_Simpson_r(1:length(Lat), 1:length(Lon)) = NaN;
Catch_Richness_p(1:length(Lat), 1:length(Lon)) = NaN;
Catch_RichnessPUE_p(1:length(Lat), 1:length(Lon)) = NaN;
Catch_Shannon_p(1:length(Lat), 1:length(Lon)) = NaN;
Catch_Simpson_p(1:length(Lat), 1:length(Lon)) = NaN;

O2_Richness_r(1:length(Lat), 1:length(Lon)) = NaN;
O2_RichnessPUE_r(1:length(Lat), 1:length(Lon)) = NaN;
O2_Shannon_r(1:length(Lat), 1:length(Lon)) = NaN;
O2_Simpson_r(1:length(Lat), 1:length(Lon)) = NaN;
O2_Richness_p(1:length(Lat), 1:length(Lon)) = NaN;
O2_RichnessPUE_p(1:length(Lat), 1:length(Lon)) = NaN;
O2_Shannon_p(1:length(Lat), 1:length(Lon)) = NaN;
O2_Simpson_p(1:length(Lat), 1:length(Lon)) = NaN;

D8_Richness_r(1:length(Lat), 1:length(Lon)) = NaN;
D8_RichnessPUE_r(1:length(Lat), 1:length(Lon)) = NaN;
D8_Shannon_r(1:length(Lat), 1:length(Lon)) = NaN;
D8_Simpson_r(1:length(Lat), 1:length(Lon)) = NaN;
D8_Richness_p(1:length(Lat), 1:length(Lon)) = NaN;
D8_RichnessPUE_p(1:length(Lat), 1:length(Lon)) = NaN;
D8_Shannon_p(1:length(Lat), 1:length(Lon)) = NaN;
D8_Simpson_p(1:length(Lat), 1:length(Lon)) = NaN;

D14_Richness_r(1:length(Lat), 1:length(Lon)) = NaN;
D14_RichnessPUE_r(1:length(Lat), 1:length(Lon)) = NaN;
D14_Shannon_r(1:length(Lat), 1:length(Lon)) = NaN;
D14_Simpson_r(1:length(Lat), 1:length(Lon)) = NaN;
D14_Richness_p(1:length(Lat), 1:length(Lon)) = NaN;
D14_RichnessPUE_p(1:length(Lat), 1:length(Lon)) = NaN;
D14_Shannon_p(1:length(Lat), 1:length(Lon)) = NaN;
D14_Simpson_p(1:length(Lat), 1:length(Lon)) = NaN;

T814_Richness_r(1:length(Lat), 1:length(Lon)) = NaN;
T814_RichnessPUE_r(1:length(Lat), 1:length(Lon)) = NaN;
T814_Shannon_r(1:length(Lat), 1:length(Lon)) = NaN;
T814_Simpson_r(1:length(Lat), 1:length(Lon)) = NaN;
T814_Richness_p(1:length(Lat), 1:length(Lon)) = NaN;
T814_RichnessPUE_p(1:length(Lat), 1:length(Lon)) = NaN;
T814_Shannon_p(1:length(Lat), 1:length(Lon)) = NaN;
T814_Simpson_p(1:length(Lat), 1:length(Lon)) = NaN;

for r = 1:1:length(Lat)
    for c = 1:1:length(Lon) 
        Catch = squeeze(Catchability(r, c, :));
        O2 = squeeze(O2_2mlpl(r, c, :));
        Rich = squeeze(Richness(r, c, :));
        RichPUE = squeeze(RichnessPUE(r, c, :));
        Shan = squeeze(Shannon(r, c, :));
        Simp = squeeze(Simpson(r, c, :));

        Depth8(1:360) = NaN;
        Depth14(1:360) = NaN;
        Thick814(1:360) = NaN;
        for m = 1:1:360
            Z8_idx = find(GODAS(r,c,:,m) >= 8);
            Z14_idx = find(GODAS(r,c,:,m) <= 14);
            
            Depth8(m) = GODAS_Depth(max(Z8_idx));
            Depth14(m) = GODAS_Depth(min(Z14_idx));
            Thick814(m) = GODAS_Depth(max(Z8_idx)) - GODAS_Depth(min(Z14_idx));
            
        end
        clear m *idx 

        % Catchability 
        [Crr, Cpr] = corrcoef(Catch, Rich, 'Rows', 'pairwise');
        Catch_Richness_r(r, c) = Crr(1,2);
        Catch_Richness_p(r, c) = Cpr(1,2);

        [Crrpue, Cprpue] = corrcoef(Catch, RichPUE, 'Rows', 'pairwise');
        Catch_RichnessPUE_r(r, c) = Crrpue(1,2);
        Catch_RichnessPUE_p(r, c) = Cprpue(1,2);

        [Crsh, Cpsh] = corrcoef(Catch, Shan, 'Rows', 'pairwise');
        Catch_Shannon_r(r, c) = Crsh(1,2);
        Catch_Shannon_p(r, c) = Cpsh(1,2);

        [Crsi, Cpsi] = corrcoef(Catch, Simp, 'Rows', 'pairwise');
        Catch_Simpson_r(r, c) = Crsi(1,2);
        Catch_Simpson_p(r, c) = Cpsi(1,2);

        % Depth of 2mlpl O2 
        [Orr, Opr] = corrcoef(O2, Rich, 'Rows', 'pairwise');
        O2_Richness_r(r, c) = Orr(1,2);
        O2_Richness_p(r, c) = Opr(1,2);

        [Orrpue, Oprpue] = corrcoef(O2, RichPUE, 'Rows', 'pairwise');
        O2_RichnessPUE_r(r, c) = Orrpue(1,2);
        O2_RichnessPUE_p(r, c) = Oprpue(1,2);
        
        [Orsh, Opsh] = corrcoef(O2, Shan, 'Rows', 'pairwise');
        O2_Shannon_r(r, c) = Orsh(1,2);
        O2_Shannon_p(r, c) = Opsh(1,2);

        [Orsi, Opsi] = corrcoef(O2, Simp, 'Rows', 'pairwise');
        O2_Simpson_r(r, c) = Orsi(1,2);
        O2_Simpson_p(r, c) = Opsi(1,2);
        
        % Depth of 8 deg C water 
        [D8rr, D8pr] = corrcoef(Depth8, Rich, 'Rows', 'pairwise');
        D8_Richness_r(r, c) = D8rr(1,2);
        D8_Richness_p(r, c) = D8pr(1,2);

        [D8rrpue, D8prpue] = corrcoef(Depth8, RichPUE, 'Rows', 'pairwise');
        D8_RichnessPUE_r(r, c) = D8rrpue(1,2);
        D8_RichnessPUE_p(r, c) = D8prpue(1,2);

        [D8rsh, D8psh] = corrcoef(Depth8, Shan, 'Rows', 'pairwise');
        D8_Shannon_r(r, c) = D8rsh(1,2);
        D8_Shannon_p(r, c) = D8psh(1,2);

        [D8rsi, D8psi] = corrcoef(Depth8, Simp, 'Rows', 'pairwise');
        D8_Simpson_r(r, c) = D8rsi(1,2);
        D8_Simpson_p(r, c) = D8psi(1,2);

        % Depth of 14 deg C water 
        [D14rr, D14pr] = corrcoef(Depth14, Rich, 'Rows', 'pairwise');
        D14_Richness_r(r, c) = D14rr(1,2);
        D14_Richness_p(r, c) = D14pr(1,2);

        [D14rrpue, D14prpue] = corrcoef(Depth14, RichPUE, 'Rows', 'pairwise');
        D14_RichnessPUE_r(r, c) = D14rrpue(1,2);
        D14_RichnessPUE_p(r, c) = D14prpue(1,2);

        [D14rsh, D14psh] = corrcoef(Depth14, Shan, 'Rows', 'pairwise');
        D14_Shannon_r(r, c) = D14rsh(1,2);
        D14_Shannon_p(r, c) = D14psh(1,2);

        [D14rsi, D14psi] = corrcoef(Depth14, Simp, 'Rows', 'pairwise');
        D14_Simpson_r(r, c) = D14rsi(1,2);
        D14_Simpson_p(r, c) = D14psi(1,2);

        % Thickness of 8-14 deg C water 
        [Trr, Tpr] = corrcoef(Thick814, Rich, 'Rows', 'pairwise');
        T814_Richness_r(r, c) = Trr(1,2);
        T814_Richness_p(r, c) = Tpr(1,2);

        [Trrpue, Tprpue] = corrcoef(Thick814, RichPUE, 'Rows', 'pairwise');
        T814_RichnessPUE_r(r, c) = Trrpue(1,2);
        T814_RichnessPUE_p(r, c) = Tprpue(1,2);

        [Trsh, Tpsh] = corrcoef(Thick814, Shan, 'Rows', 'pairwise');
        T814_Shannon_r(r, c) = Trsh(1,2);
        T814_Shannon_p(r, c) = Tpsh(1,2);

        [Trsi, Tpsi] = corrcoef(Thick814, Simp, 'Rows', 'pairwise');
        T814_Simpson_r(r, c) = Trsi(1,2);
        T814_Simpson_p(r, c) = Tpsi(1,2);

    end
end
clear r c CP* DP* OP* Depth8 Depth14 Thick814 Catch O2 Hooks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plots
% There are functions for this at the end of the script

% Shannon Index
figure
Corr_Map(lat_grid, lon_grid, Catch_Shannon_r, Catch_Shannon_p, Vessels_total, ...
    0.05, -1, 1, "Pearson Correlation: Catchability and Shannon Index")

figure
Corr_Map(lat_grid, lon_grid, O2_Shannon_r, O2_Shannon_p, Vessels_total, ...
    0.05, -1, 1, "Pearson Correlation: Depth of 2ml l^-^1 O_2 Surface and Shannon Index")

figure
Corr_Map(lat_grid, lon_grid, D8_Shannon_r, D8_Shannon_p, Vessels_total, ...
    0.05, -1, 1, "Pearson Correlation: Depth of 8-deg Surface and Shannon Index")

figure
Corr_Map(lat_grid, lon_grid, D14_Shannon_r, D14_Shannon_p, Vessels_total, ...
    0.05, -1, 1, "Pearson Correlation: Depth of 14-deg Surface and Shannon Index")

figure
Corr_Map(lat_grid, lon_grid, T814_Shannon_r, T814_Shannon_p, Vessels_total, ...
    0.05, -1, 1, "Pearson Correlation: Thickness of 8–14-deg Waters and Shannon Index")

% Simpson Index
figure
Corr_Map(lat_grid, lon_grid, Catch_Simpson_r, Catch_Simpson_p, Vessels_total, ...
    0.05, -1, 1, "Pearson Correlation: Catchability and Simpson Index")

figure
Corr_Map(lat_grid, lon_grid, O2_Simpson_r, O2_Simpson_p, Vessels_total, ...
    0.05, -1, 1, "Pearson Correlation: Depth of 2ml l^-^1 O_2 Surface and Simpson Index")

figure
Corr_Map(lat_grid, lon_grid, D8_Simpson_r, D8_Simpson_p, Vessels_total, ...
    0.05, -1, 1, "Pearson Correlation: Depth of 8-deg Surface and Simpson Index")

figure
Corr_Map(lat_grid, lon_grid, D14_Simpson_r, D14_Simpson_p, Vessels_total, ...
    0.05, -1, 1, "Pearson Correlation: Depth of 14-deg Surface and Simpson Index")

figure
Corr_Map(lat_grid, lon_grid, T814_Simpson_r, T814_Simpson_p, Vessels_total, ...
    0.05, -1, 1, "Pearson Correlation: Thickness of 8–14-deg Waters and Simpson Index")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Functions
% Function to create map of significant correlations and metrics themselves
function Corr_Map(lat_grid_used, lon_grid_used, r_vals, p_vals, vessels, sig_level, min_val, max_val, map_title)
r_vals(p_vals >= sig_level) = NaN;
r_vals(vessels < 3) = NaN;
axesm('mercator','MapLatLimit',[10 40],'MapLonLimit',[180 230], ...
    'MLineLocation', 10, 'PLineLocation', 10, ... % draw every 10 degrees         
    'Grid', 'on', 'MeridianLabel','on','ParallelLabel','on', ...
    'MLabelParallel', 10, 'MLabelLocation', 10, 'PLabelLocation', 10); % label every 10 degrees, below map 
Ig = geoshow(lat_grid_used,lon_grid_used,r_vals,'DisplayType','texturemap');
set(Ig,'AlphaData',double(~isnan(r_vals)),'AlphaDataMapping','none','FaceAlpha','texturemap'); 
plotm([20 20], [180 230], 'k');
plotm([10 40], [210 210], 'k');
plotm([26 26], [180 210], 'k');
geoshow('landareas.shp','FaceColor',[0.5 0.5 0.5]);
colormap(redblueTecplot)
clim([min_val max_val])
cb = colorbar();
ylabel(cb, 'Correlation Coefficient', 'Rotation', 270);
title(sprintf('%s'), map_title)
set(gcf,'renderer','Painters')
tightmap
end






























