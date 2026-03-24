% The purpose of this script is to examine whether there's a significant
% linear trend in effort in each grid cell.

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
% PC:
% Lat = ncread('..\FisheryData\TotalEffort.nc', 'Latitude');
% Lon = ncread('..\FisheryData\TotalEffort.nc', 'Longitude');
% Month = ncread('..\FisheryData\TotalEffort.nc', 'Month');
% Effort = ncread('..\FisheryData\TotalEffort.nc', 'Total Effort');
% Vessels = ncread('..\FisheryData\TotalVessels.nc', 'Total Vessels');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Wrangle data
% Sum Effort each of the given regions
% Permute so arrays are Lat x Lon x time
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
%%% Linear trends
% Predictor: sequential months, 1–360
x = 1:1:360;

% Empty matrix to fill
EffortTrends(1:length(Lat), 1:length(Lon)) = NaN;

% Loop through all grid cells
for r = 1:1:length(Lat)
    for c = 1:1:length(Lon)
        
        % Only model cells with at least 5 years of effort (arbitrary)
        if sum(Effort(r,c,:) > 0) > 60

            % Linear model
            px_lm = fitlm(x, squeeze(Effort(r,c,:)));

            % Assess p value and fill map with slope if < 0.05
            if px_lm.Coefficients.pValue(2) < 0.05

                % Calculating the slope to avoid wrangling a table
                EffortTrends(r,c) = round((px_lm.Fitted(end) - px_lm.Fitted(1)) / (x(end) - x(1)));
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Map
figure
axesm('mercator','MapLatLimit',[10 40],'MapLonLimit',[180 230], ...
    'MLineLocation', 10, 'PLineLocation', 10, ... % draw every 10 degrees         
    'Grid', 'on', 'MeridianLabel','on','ParallelLabel','on', ...
    'MLabelParallel', 10, 'MLabelLocation', 10, 'PLabelLocation', 10); % label every 10 degrees, below map 
Ig = geoshow(lat_grid,lon_grid,EffortTrends,'DisplayType','texturemap');
set(Ig,'AlphaData',double(~isnan(EffortTrends)),'AlphaDataMapping','none','FaceAlpha','texturemap'); 
% plotm([20 20], [180 230], 'k');
% plotm([10 40], [210 210], 'k');
% plotm([26 26], [180 210], 'k');
geoshow('landareas.shp','FaceColor',[0.5 0.5 0.5]);
colormap(redblueTecplot)
clim([-180 180])
cb = colorbar();
ylabel(cb, 'Hooks per Month', 'Rotation', 270);
title('Significant Linear Effort Trends')
set(gcf,'renderer','Painters')
tightmap
