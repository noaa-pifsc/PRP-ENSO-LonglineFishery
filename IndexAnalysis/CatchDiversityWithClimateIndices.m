% The purpose of this script is to look at how measures of catch diversity 
% time in each of the regions of the fishing grounds we're examining and to
% compare this with the various climate indices we're using. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Load data
% Matlab needs a "here" package.  Until then, use the path that fits your
% OS.
% Mac:
Lat = ncread('../FisheryData/TotalEffort.nc', 'Latitude');
Lon = ncread('../FisheryData/TotalEffort.nc', 'Longitude');
Month = ncread('../FisheryData/TotalEffort.nc', 'Month');
Effort = ncread('../FisheryData/TotalEffort.nc', 'Total Effort');
Vessels = ncread('../FisheryData/TotalVessels.nc', 'Total Vessels');
ONI = readtable('../ClimateIndices/ONI_withPhases.csv');
PDO = readtable('../ClimateIndices/PDO.csv');
NPGO = readtable('../ClimateIndices/NPGO.csv');
Richness = ncread('../FisheryData/SpeciesRichness.nc', 'Species Richness');
Shannon = ncread('../FisheryData/ShannonIndex.nc', 'Shannon Index');
Simpson = ncread('../FisheryData/SimpsonIndex.nc', 'Simpson Index');
% PC:
% Lat = ncread('..\FisheryData\TotalEffort.nc', 'Latitude');
% Lon = ncread('..\FisheryData\TotalEffort.nc', 'Longitude');
% Month = ncread('..\FisheryData\TotalEffort.nc', 'Month');
% Effort = ncread('..\FisheryData\TotalEffort.nc', 'Total Effort');
% Vessels = ncread('..\FisheryData\TotalVessels.nc', 'Total Vessels');
% ONI = readtable('..\ClimateIndices\ONI_withPhases.csv');
% PDO = readtable('..\ClimateIndices\PDO.csv');
% NPGO = readtable('..\ClimateIndices\NPGO.csv');
% Richness = ncread('..\FisheryData\SpeciesRichness.nc', 'Species Richness');
% Shannon = ncread('..\FisheryData\ShannonIndex.nc', 'Shannon Index');
% Simpson = ncread('..\FisheryData\SimpsonIndex.nc', 'Simpson Index');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Wrangle data
% Limit climate indices to our period of interest: 1995–2024
ONI_9524 = ONI(ONI.YR >= 1995 & ONI.YR <= 2024,:);
PDO_9524 = PDO(PDO.Year >= 1995 & PDO.Year <= 2024,:);
NPGO_9524 = NPGO(NPGO.YEAR >= 1995 & NPGO.YEAR <= 2024,:);
% Clean up
clear ONI PDO NPGO

% To make nice plots, it is helpful to separate out positive and negative
% phases of each index (which is a little tedious)
PosONI = ONI_9524.ONI;
PosONI(isnan(str2double(ONI_9524.ElNino))) = NaN;
NegONI = ONI_9524.ONI;
NegONI(isnan(str2double(ONI_9524.LaNina))) = NaN;
NeuONI = ONI_9524.ONI;
NeuONI(isnan(str2double(ONI_9524.Neutral))) = NaN;
NonONI = ONI_9524.ONI;
NonONI(~isnan(str2double(ONI_9524.ElNino)) | ~isnan(str2double(ONI_9524.LaNina)) | ~isnan(str2double(ONI_9524.Neutral))) = NaN;
ONI_9524 = addvars(ONI_9524, PosONI);
ONI_9524 = addvars(ONI_9524, NegONI);
ONI_9524 = addvars(ONI_9524, NeuONI);
ONI_9524 = addvars(ONI_9524, NonONI);

PosPDO = PDO_9524.PDO;
PosPDO(PosPDO < 0) = NaN;
NegPDO = PDO_9524.PDO;
NegPDO(NegPDO > 0) = NaN;
NeuPDO = PDO_9524.PDO;
NeuPDO(NeuPDO ~= 0) = NaN;
PDO_9524 = addvars(PDO_9524, PosPDO);
PDO_9524 = addvars(PDO_9524, NegPDO);
PDO_9524 = addvars(PDO_9524, NeuPDO);

PosNPGO = NPGO_9524.NPGO;
PosNPGO(PosNPGO < 0) = NaN;
NegNPGO = NPGO_9524.NPGO;
NegNPGO(NegNPGO > 0) = NaN;
NeuNPGO = NPGO_9524.NPGO;
NeuNPGO(NeuNPGO ~= 0) = NaN;
NPGO_9524 = addvars(NPGO_9524, PosNPGO);
NPGO_9524 = addvars(NPGO_9524, NegNPGO);

% Permute so arrays are Lat x Lon x time
Effort = permute(Effort, [2 1 3]);
Richness = permute(Richness, [2 1 3]);
Shannon = permute(Shannon, [2 1 3]);
Simpson = permute(Simpson, [2 1 3]);
Vessels = permute(Vessels, [2 1 3]);

% Average over space each year
RichnessTS = squeeze(mean(Richness, [1 2], "omitnan"));
ShannonTS = squeeze(mean(Shannon, [1 2], "omitnan"));
SimpsonTS = squeeze(mean(Simpson, [1 2], "omitnan"));

% Sum over domain each year
EffortTS = squeeze(sum(Effort, [1 2], "omitnan"));

% Climatologies/totals
Richness_climo = mean(Richness, 3, "omitnan");
Shannon_climo = mean(Shannon, 3, "omitnan");
Simpson_climo = mean(Simpson, 3, "omitnan");
Effort_total = sum(Effort, 3, "omitnan");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Correlations
% Compare to climate variability index
[O_rich_r, O_rich_p] = corrcoef(RichnessTS(:,1), ONI_9524.ONI, 'Rows', 'pairwise');
[O_shan_r, O_shan_p] = corrcoef(ShannonTS(:,1), ONI_9524.ONI, 'Rows', 'pairwise');
[O_simp_r, O_simp_p] = corrcoef(SimpsonTS(:,1), ONI_9524.ONI, 'Rows', 'pairwise');

[P_rich_r, P_rich_p] = corrcoef(RichnessTS(:,1), PDO_9524.PDO, 'Rows', 'pairwise');
[P_shan_r, P_shan_p] = corrcoef(ShannonTS(:,1), PDO_9524.PDO, 'Rows', 'pairwise');
[P_simp_r, P_simp_p] = corrcoef(SimpsonTS(:,1), PDO_9524.PDO, 'Rows', 'pairwise');

[N_rich_r, N_rich_p] = corrcoef(RichnessTS(:,1), NPGO_9524.NPGO, 'Rows', 'pairwise');
[N_shan_r, N_shan_p] = corrcoef(ShannonTS(:,1), NPGO_9524.NPGO, 'Rows', 'pairwise');
[N_simp_r, N_simp_p] = corrcoef(SimpsonTS(:,1), NPGO_9524.NPGO, 'Rows', 'pairwise');

% Account for relationship between richness and effort
RPUE = RichnessTS ./ EffortTS * 10^6;
[O_richpue_r, O_richpue_p] = corrcoef(RPUE(:,1), ONI_9524.ONI, 'Rows', 'pairwise');
[P_richpue_r, P_richpue_p] = corrcoef(RPUE(:,1), PDO_9524.PDO, 'Rows', 'pairwise');
[N_richpue_r, N_richpue_p] = corrcoef(RPUE(:,1), NPGO_9524.NPGO, 'Rows', 'pairwise');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Time series plots for significant correlations
% There are functions for this at the end of the script
% PDO_yyPlot(climate, div_metric, metric_name)
% NPGO_yyPlot(climate, div_metric, metric_name)

% figure
% PDO_yyPlot(PDO_9524, RichnessTS, 'Species Richness');
% 
% figure
% PDO_yyPlot(PDO_9524, RichnessTS./EffortTS*10^6, 'Species Richness per Unit Effort');
% 
% figure
% PDO_yyPlot(PDO_9524, ShannonTS, 'Shannon Index');
% 
% figure
% PDO_yyPlot(PDO_9524, SimpsonTS, 'Simpson Index');

% figure
% NPGO_yyPlot(NPGO_9524, RichnessTS, 'Species Richness');
% 
% figure
% NPGO_yyPlot(NPGO_9524, RichnessTS./EffortTS*10^6, 'Species Richness per Unit Effort');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Create maps of richness/diversity by phase
% Create lat and lon grids
% Horizontally stack column with latitudes
lat_grid = repmat(Lat,1,size(Lon,1));
% Transpose column with longitudes to a row and vertically stack
lon_grid = repmat(Lon',size(Lat,1),1);

% Create a mask to make maps non-confidential
Tot_vessels = sum(Vessels, 3, 'omitnan');
Confid_mask(1:size(Vessels, 1), 1:size(Vessels, 2)) = ones;
Confid_mask(Tot_vessels < 3) = NaN;

% There are functions for this at the end of the script
% Metric_Map(metric_to_map, map_mask, lat_grid_used, lon_grid_used, climate_mode_phase, mn, mx, panel_title)
% Climo_Map(metric_to_map, map_mask, lat_grid_used, lon_grid_used, panel_title)
% Climatologies
% figure
% Climo_Map(Richness_climo, Confid_mask, lat_grid, lon_grid, "Species Richness")
% 
% figure
% Climo_Map(Richness_climo./Effort_total*10^6, Confid_mask, lat_grid, lon_grid, "Species Richness per Unit Effort")

% figure
% Climo_Map(Shannon_climo, Confid_mask, lat_grid, lon_grid, "Shannon Index")
% 
% figure
% Climo_Map(Simpson_climo, Confid_mask, lat_grid, lon_grid, "Simpson Index")

% % ONI
% figure
% subplot(2,2,1)
% Metric_Map(Richness, Confid_mask, lat_grid, lon_grid, PosONI, 1, 15, 'El Niño')
% 
% subplot(2,2,2)
% Metric_Map(Richness, Confid_mask, lat_grid, lon_grid, NegONI, 1, 15, 'La Niña')
% 
% subplot(2,2,3)
% Metric_Map(Richness, Confid_mask, lat_grid, lon_grid, NeuONI, 1, 15, 'Neutral')
% 
% subplot(2,2,4)
% Metric_Map(Richness, Confid_mask, lat_grid, lon_grid, NonONI, 1, 15, '<5 Consecutive Months')
% 
% PDO
figure
subplot(1,2,1)
Metric_Map(Shannon, Confid_mask, lat_grid, lon_grid, PosPDO, 0, 2.8, 'Shannon Index in Positive PDO')

subplot(1,2,2)
Metric_Map(Shannon, Confid_mask, lat_grid, lon_grid, NegPDO, 0, 2.8, 'Shannon Index in Negative PDO')


figure
subplot(1,2,1)
Metric_Map(Simpson, Confid_mask, lat_grid, lon_grid, PosPDO, 0, 1, 'Simpson Index in Positive PDO')

subplot(1,2,2)
Metric_Map(Simpson, Confid_mask, lat_grid, lon_grid, NegPDO, 0, 1, 'Simpson Index in Negative PDO')

% % NPGO
% figure
% subplot(1,2,1)
% Metric_Map(Shannon, Confid_mask, lat_grid, lon_grid, PosNPGO, 1, 15, 'Positive NPGO')
% 
% subplot(1,2,2)
% Metric_Map(Shannon, Confid_mask, lat_grid, lon_grid, NegNPGO, 1, 15, 'Negative NPGO')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Functions
% Function to create richness + PDO plots
function PDO_yyPlot(climate, div_metric, metric_name)
colororder({'#A5AAAF', 'k'}) % y-axes colors, left to right
yyaxis left
stem(climate.PosPDO, '-r', 'Marker', 'none');
hold on
stem(climate.NegPDO, '-b', 'Marker', 'none');
ylim([-3 3]);
ylabel('Pacific Decadal Oscillation');
yyaxis right
plot(div_metric);
ylabel(sprintf('%s', metric_name));
xlim([1 360]);
set(gca,'XTick',1:60:360);
set(gca,'XTickLabel',1995:5:2020);
pbaspect([2 1 1]);
end


% Function to create effort + NPGO plots
function NPGO_yyPlot(climate, div_metric, metric_name)
colororder({'#A5AAAF', 'k'}) % y-axes colors, left to right
yyaxis left
stem(climate.PosNPGO, '-r', 'Marker', 'none');
hold on
stem(climate.NegNPGO, '-b', 'Marker', 'none');
ylim([-3.5 3.5]);
ylabel('North Pacific Gyre Oscillation');
yyaxis right
plot(div_metric);
ylabel(sprintf('%s', metric_name));
xlim([1 360]);
set(gca,'XTick',1:60:360);
set(gca,'XTickLabel',1995:5:2020);
pbaspect([2 1 1]);
end

% Function to create maps of diversity metrics by phase
function Metric_Map(metric_to_map, map_mask, lat_grid_used, lon_grid_used, climate_mode_phase, mn, mx, panel_title)
m_map = mean(metric_to_map(:,:,~isnan(climate_mode_phase)),3, "omitnan");
m_map(m_map == 0) = NaN; % zeros clutter the plot
m_map = m_map .* map_mask; % to eliminate confidential points
axesm('mercator','MapLatLimit',[10 40],'MapLonLimit',[180 230], ...
    'MLineLocation', 10, 'PLineLocation', 10, ... % draw every 10 degrees         
    'Grid', 'on', 'MeridianLabel','on','ParallelLabel','on', ...
    'MLabelParallel', 10, 'MLabelLocation', 10, 'PLabelLocation', 10); % label every 10 degrees, below map 
Ig = geoshow(lat_grid_used,lon_grid_used,m_map,'DisplayType','texturemap');
set(Ig,'AlphaData',double(~isnan(m_map)),'AlphaDataMapping','none','FaceAlpha','texturemap'); 
plotm([20 20], [180 230], 'k');
plotm([10 40], [210 210], 'k');
plotm([26 26], [180 210], 'k');
geoshow('landareas.shp','FaceColor',[0.5 0.5 0.5]);
clim([mn mx])
cb = colorbar();
ylabel(cb, ' ', 'Rotation', 270);
title(sprintf('%s'), panel_title)
set(gcf,'renderer','Painters')
tightmap
end

% Function to create climatological maps of diversity metrics
function Climo_Map(metric_to_map, map_mask, lat_grid_used, lon_grid_used, panel_title)
metric_to_map(metric_to_map == 0) = NaN; % zeros clutter the map
metric_to_map = metric_to_map .* map_mask; % to eliminate confidential points
axesm('mercator','MapLatLimit',[10 40],'MapLonLimit',[180 230], ...
    'MLineLocation', 10, 'PLineLocation', 10, ... % draw every 10 degrees         
    'Grid', 'on', 'MeridianLabel','on','ParallelLabel','on', ...
    'MLabelParallel', 10, 'MLabelLocation', 10, 'PLabelLocation', 10); % label every 10 degrees, below map 
Ig = geoshow(lat_grid_used,lon_grid_used,metric_to_map,'DisplayType','texturemap');
set(Ig,'AlphaData',double(~isnan(metric_to_map)),'AlphaDataMapping','none','FaceAlpha','texturemap'); 
plotm([20 20], [180 230], 'k');
plotm([10 40], [210 210], 'k');
plotm([26 26], [180 210], 'k');
geoshow('landareas.shp','FaceColor',[0.5 0.5 0.5]);
cb = colorbar();
ylabel(cb, ' ', 'Rotation', 270);
title(sprintf('%s'), panel_title)
set(gcf,'renderer','Painters')
tightmap
end












