% The purpose of this script is to look at how effort is distributed over
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
% PC:
% Lat = ncread('..\FisheryData\TotalEffort.nc', 'Latitude');
% Lon = ncread('..\FisheryData\TotalEffort.nc', 'Longitude');
% Month = ncread('..\FisheryData\TotalEffort.nc', 'Month');
% Effort = ncread('..\FisheryData\TotalEffort.nc', 'Total Effort');
% Vessels = ncread('..\FisheryData\TotalVessels.nc', 'Total Vessels');
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


% Sum Effort each of the given regions
% Permute so arrays are Lat x Lon x time
Effort = permute(Effort, [2 1 3]);
TotalEffort = squeeze(sum(Effort, [1 2], "omitnan"));

% Separate Effort into regions
NW_lat = find(Lat >= 26);
CW_lat = find(Lat >= 20 & Lat < 26);
S_lat = find(Lat >= 10 & Lat < 20);
NE_lat = find(Lat >= 20);
W_lon = find(Lon >= 180 & Lon < 210);
E_lon = find (Lon >= 210 & Lon < 230);

Regional_Effort(1:length(Month), 1:5) = NaN; %Cols are NW, CW, SW, SE, NE
for m = 1:1:length(Month)
    Regional_Effort(m,1) = sum(Effort(NW_lat, W_lon, m), "all", "omitnan");
    Regional_Effort(m,2) = sum(Effort(CW_lat, W_lon, m), "all", "omitnan");
    Regional_Effort(m,3) = sum(Effort(S_lat, W_lon, m), "all", "omitnan");
    Regional_Effort(m,4) = sum(Effort(S_lat, E_lon, m), "all", "omitnan");
    Regional_Effort(m,5) = sum(Effort(NE_lat, E_lon, m), "all", "omitnan");
end
clear m NW_lat CW_lat S_lat NE_lat W_lon E_lon

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Correlations
% Compare to climate variability index
O_Spearman_r(1,1:5) = NaN;
O_Spearman_p(1,1:5) = NaN;
P_Spearman_r(1,1:5) = NaN;
P_Spearman_p(1,1:5) = NaN;
N_Spearman_r(1,1:5) = NaN;
N_Spearman_p(1,1:5) = NaN;

O_Pearson_r(1,1:5) = NaN;
O_Pearson_p(1,1:5) = NaN;
P_Pearson_r(1,1:5) = NaN;
P_Pearson_p(1,1:5) = NaN;
N_Pearson_r(1,1:5) = NaN;
N_Pearson_p(1,1:5) = NaN;

for r = 1:1:5
    [O_Spearman_r(1,r), O_Spearman_p(1,r)] = corr(Regional_Effort(:,r), ONI_9524.ONI,...
        'type', 'Spearman', 'rows', 'complete');
    [P_Spearman_r(1,r), P_Spearman_p(1,r)] = corr(Regional_Effort(:,r), PDO_9524.PDO,...
        'type', 'Spearman', 'rows', 'complete');
    [N_Spearman_r(1,r), N_Spearman_p(1,r)] = corr(Regional_Effort(:,r), NPGO_9524.NPGO,...
        'type', 'Spearman', 'rows', 'complete');

    [Or, Op] = corrcoef(Regional_Effort(:,r), ONI_9524.ONI, 'Rows', 'pairwise');
    [Pr, Pp] = corrcoef(Regional_Effort(:,r), PDO_9524.PDO, 'Rows', 'pairwise');
    [Nr, Np] = corrcoef(Regional_Effort(:,r), NPGO_9524.NPGO, 'Rows', 'pairwise');

    O_Pearson_r(1,r) = Or(1,2);
    O_Pearson_p(1,r) = Op(1,2);
    P_Pearson_r(1,r) = Pr(1,2);
    P_Pearson_p(1,r) = Pp(1,2);
    N_Pearson_r(1,r) = Nr(1,2);
    N_Pearson_p(1,r) = Np(1,2);
end
clear r

[O_tot_r, O_tot_p] = corrcoef(TotalEffort(:,1), ONI_9524.ONI, 'Rows', 'pairwise');
[P_tot_r, P_tot_p] = corrcoef(TotalEffort(:,1), PDO_9524.PDO, 'Rows', 'pairwise');
[N_tot_r, N_tot_p] = corrcoef(TotalEffort(:,1), NPGO_9524.NPGO, 'Rows', 'pairwise');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Time series plots for significant correlations
% There are functions for this at the end of the script
% ONI_Effort_yyPlot(region_name, region_number, climate, fishery)
% PDO_Effort_yyPlot(region_name, region_number, climate, fishery)
% NPGO_Effort_yyPlot(region_name, region_number, climate, fishery)
figure
ONI_Effort_yyPlot('Centralwest', 2, ONI_9524, Regional_Effort);

figure
PDO_Effort_yyPlot('Northwest', 1, PDO_9524, Regional_Effort);

figure
PDO_Effort_yyPlot('Centralwest', 2, PDO_9524, Regional_Effort);

figure
PDO_Effort_yyPlot('Southwest', 3, PDO_9524, Regional_Effort);

figure
PDO_Effort_yyPlot('Southeast', 4, PDO_9524, Regional_Effort);

figure
PDO_Effort_yyPlot('Northeast', 5, PDO_9524, Regional_Effort);

figure
NPGO_Effort_yyPlot('Northwest', 1, NPGO_9524, Regional_Effort);

figure
NPGO_Effort_yyPlot('Centralwest', 2, NPGO_9524, Regional_Effort);

figure
NPGO_Effort_yyPlot('Southwest', 3, NPGO_9524, Regional_Effort);

figure
NPGO_Effort_yyPlot('Southeast', 4, NPGO_9524, Regional_Effort);

figure
NPGO_Effort_yyPlot('Northeast', 5, NPGO_9524, Regional_Effort);

% Total Effort
figure
colororder({'#A5AAAF', 'k'}) % y-axes colors, left to right
yyaxis left
stem(PDO_9524.PosPDO, '-r', 'Marker', 'none');
hold on
stem(PDO_9524.NegPDO, '-b', 'Marker', 'none');
ylim([-3 3]);
ylabel('Pacific Decadal Oscillation');
yyaxis right
plot(TotalEffort(:,1));
ylim([0  6e6]);
ylabel('Hooks');
xlim([1 360]);
set(gca,'XTick',1:60:360);
set(gca,'XTickLabel',1995:5:2020);
title('Full Fishery')
pbaspect([2 1 1]);

figure
colororder({'#A5AAAF', 'k'}) % y-axes colors, left to right
yyaxis left
stem(NPGO_9524.PosNPGO, '-r', 'Marker', 'none');
hold on
stem(NPGO_9524.NegNPGO, '-b', 'Marker', 'none');
ylim([-3.5 3.5]);
ylabel('North Pacific Gyre Oscillation');
yyaxis right
plot(TotalEffort(:,1));
ylim([0  6e6]);
ylabel('Hooks');
xlim([1 360]);
set(gca,'XTick',1:60:360);
set(gca,'XTickLabel',1995:5:2020);
title('Full Fishery')
pbaspect([2 1 1]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Comparison of how months and effort are distributed across phases
% Create arrays for each mode that are month x phase x region (6 = full
% fishery)
ONI_array(1:size(ONI_9524,1),1:4,1:6) = NaN; % pos, neg, neutral, non
PDO_array(1:size(PDO_9524,1),1:2,1:6) = NaN; % pos, neg
NPGO_array(1:size(NPGO_9524,1),1:2,1:6) = NaN; % pos, neg

% Add total effort to regional effort to make plotting easier
Regional_Effort(:,6) = TotalEffort(:,1);

% Loop through regions to fill
for r = 1:1:6
    % Find months with effort in given region
    mo_idx = find(Regional_Effort(:,r) > 0);
    
    % Fill arrays with index for months that have effort
    % ONI
    ONI_array(mo_idx,1,r) = PosONI(mo_idx,1);
    ONI_array(mo_idx,2,r) = NegONI(mo_idx,1);
    ONI_array(mo_idx,3,r) = NeuONI(mo_idx,1);
    ONI_array(mo_idx,4,r) = NonONI(mo_idx,1);

    % PDO
    PDO_array(mo_idx,1,r) = PosPDO(mo_idx,1);
    PDO_array(mo_idx,2,r) = NegPDO(mo_idx,1);
    
    % NPGO
    NPGO_array(mo_idx,1,r) = PosNPGO(mo_idx,1);
    NPGO_array(mo_idx,2,r) = NegNPGO(mo_idx,1); 
end
clear r mo_idx

% Grids for plotting
ONI_Grid(1:9,1) = 10:10:90;
ONI_Grid(1:9,2) = 10:10:90;
ONI_Grid(1:9,3) = 0.2;
ONI_Grid(1:9,4) = 6.6;

Other_Grid(1:9,1) = 10:10:90;
Other_Grid(1:9,2) = 10:10:90;
Other_Grid(1:9,3) = 0.6;
Other_Grid(1:9,4) = 6.6;

% Plot
figure
for r = 1:1:6
    % Thick line for percent of months with effort in a given phase 
    plot([0 sum(~isnan(ONI_array(:,1,r))) / sum(Regional_Effort(:,r) > 0) * 100],...
        [r r], 'LineWidth', 7, 'Color', '#FFBFBF'); % El Nino
    hold on
    plot([0 sum(~isnan(ONI_array(:,2,r))) / sum(Regional_Effort(:,r) > 0) * 100],...
        [r+0.2 r+0.2], 'LineWidth', 7, 'Color', '#BFBFFF'); % La Nina
    hold on
    plot([0 sum(~isnan(ONI_array(:,3,r))) / sum(Regional_Effort(:,r) > 0) * 100],...
        [r+0.4 r+0.4], 'LineWidth', 7, 'Color', '#BFBFBF'); % Neutral
    hold on
    plot([0 sum(~isnan(ONI_array(:,4,r))) / sum(Regional_Effort(:,r) > 0) * 100],...
        [r+0.6 r+0.6], 'LineWidth', 7, 'Color', '#DFDFDF'); % None
    hold on

    % Markers for the percent of months in a given phase with effort
    plot(sum(~isnan(ONI_array(:,1,r))) / sum(~isnan(PosONI(:,1))) * 100,...
        r, 'o', 'MarkerFaceColor', '#FFBFBF', 'MarkerEdgeColor', '#FF1e1e');
    hold on
    plot(sum(~isnan(ONI_array(:,2,r))) / sum(~isnan(NegONI(:,1))) * 100,...
        r+0.2, 'o', 'MarkerFaceColor', '#BFBFFF', 'MarkerEdgeColor', '#1111FF');
    hold on
    plot(sum(~isnan(ONI_array(:,3,r))) / sum(~isnan(NeuONI(:,1))) * 100,...
        r+0.4, 'o', 'MarkerFaceColor', '#BFBFBF', 'MarkerEdgeColor', '#151515');
    hold on
    plot(sum(~isnan(ONI_array(:,4,r))) / sum(~isnan(NonONI(:,1))) * 100,...
        r+0.6, 'o', 'MarkerFaceColor', '#DFDFDF', 'MarkerEdgeColor', '#929292');
    hold on

    % Thin dark line for percent of effort set in a given region in a given phase 
    plot([0 sum(Regional_Effort(~isnan(ONI_array(:,1,r)),r)) / sum(Regional_Effort(:,r)) * 100],...
        [r r], 'Color', '#FF1e1e');
    hold on
    plot([0 sum(Regional_Effort(~isnan(ONI_array(:,2,r)),r)) / sum(Regional_Effort(:,r)) * 100],...
        [r+0.2 r+0.2], 'Color', '#1111FF');
    hold on
    plot([0 sum(Regional_Effort(~isnan(ONI_array(:,3,r)),r)) / sum(Regional_Effort(:,r)) * 100],...
        [r+0.4 r+0.4], 'Color', '#151515');
    hold on
    plot([0 sum(Regional_Effort(~isnan(ONI_array(:,4,r)),r)) / sum(Regional_Effort(:,r)) * 100],...
        [r+0.6 r+0.6], 'Color', '#929292');
    hold on
end
% Grid
for l = 1:1:9
    plot([ONI_Grid(l,1) ONI_Grid(l,2)], [ONI_Grid(l,3) ONI_Grid(l,4)], 'w');
    hold on
end
set(gca, 'YDir', 'Reverse'); 
set(gca, 'YTick', 1.3:1:6.3);
set(gca, 'YTickLabel', {'NW', 'CW', 'SW', 'SE', 'NE', 'All'});
xlabel('Percent');
axis([0 100 0.6 7]);
title('Oceanic Niño Index');
pbaspect([2 1 1]);
clear r l

figure
for r = 1:1:6
    % Thick line for percent of months with effort in a given phase 
    plot([0 sum(~isnan(PDO_array(:,1,r))) / sum(Regional_Effort(:,r) > 0) * 100],...
        [r r], 'LineWidth', 7, 'Color', '#FFBFBF'); % Positive
    hold on
    plot([0 sum(~isnan(PDO_array(:,2,r))) / sum(Regional_Effort(:,r) > 0) * 100],...
        [r+0.2 r+0.2], 'LineWidth', 7, 'Color', '#BFBFFF'); % Negative
    hold on
    
    % Markers for the percent of months in a given phase with effort
    plot(sum(~isnan(PDO_array(:,1,r))) / sum(~isnan(PosPDO(:,1))) * 100,...
        r, 'o', 'MarkerFaceColor', '#FFBFBF', 'MarkerEdgeColor', '#FF1e1e');
    hold on
    plot(sum(~isnan(PDO_array(:,2,r))) / sum(~isnan(NegPDO(:,1))) * 100,...
        r+0.2, 'o', 'MarkerFaceColor', '#BFBFFF', 'MarkerEdgeColor', '#1111FF');
    hold on
    
    % Thin dark line for percent of effort set in a given region in a given phase 
    plot([0 sum(Regional_Effort(~isnan(PDO_array(:,1,r)),r)) / sum(Regional_Effort(:,r)) * 100],...
        [r r], 'Color', '#FF1e1e');
    hold on
    plot([0 sum(Regional_Effort(~isnan(PDO_array(:,2,r)),r)) / sum(Regional_Effort(:,r)) * 100],...
        [r+0.2 r+0.2], 'Color', '#1111FF');
    hold on
end
% Grid
for l = 1:1:9
    plot([Other_Grid(l,1) Other_Grid(l,2)], [Other_Grid(l,3) Other_Grid(l,4)], 'w');
    hold on
end
set(gca, 'YDir', 'Reverse'); 
set(gca, 'YTick', 1.1:1:6.1);
set(gca, 'YTickLabel', {'NW', 'CW', 'SW', 'SE', 'NE', 'All'});
xlabel('Percent');
axis([0 100 0.2 7]);
title('Pacific Decadal Index');
pbaspect([2 1 1]);
clear r l

figure
for r = 1:1:6
    % Thick line for percent of months with effort in a given phase 
    plot([0 sum(~isnan(NPGO_array(:,1,r))) / sum(Regional_Effort(:,r) > 0) * 100],...
        [r r], 'LineWidth', 7, 'Color', '#FFBFBF'); % Positive
    hold on
    plot([0 sum(~isnan(NPGO_array(:,2,r))) / sum(Regional_Effort(:,r) > 0) * 100],...
        [r+0.2 r+0.2], 'LineWidth', 7, 'Color', '#BFBFFF'); % Negative
    hold on
    
    % Markers for the percent of months in a given phase with effort
    plot(sum(~isnan(NPGO_array(:,1,r))) / sum(~isnan(PosNPGO(:,1))) * 100,...
        r, 'o', 'MarkerFaceColor', '#FFBFBF', 'MarkerEdgeColor', '#FF1e1e');
    hold on
    plot(sum(~isnan(NPGO_array(:,2,r))) / sum(~isnan(NegNPGO(:,1))) * 100,...
        r+0.2, 'o', 'MarkerFaceColor', '#BFBFFF', 'MarkerEdgeColor', '#1111FF');
    hold on
    
    % Thin dark line for percent of effort set in a given region in a given phase 
    plot([0 sum(Regional_Effort(~isnan(NPGO_array(:,1,r)),r)) / sum(Regional_Effort(:,r)) * 100],...
        [r r], 'Color', '#FF1e1e');
    hold on
    plot([0 sum(Regional_Effort(~isnan(NPGO_array(:,2,r)),r)) / sum(Regional_Effort(:,r)) * 100],...
        [r+0.2 r+0.2], 'Color', '#1111FF');
    hold on
end
% Grid
for l = 1:1:9
    plot([Other_Grid(l,1) Other_Grid(l,2)], [Other_Grid(l,3) Other_Grid(l,4)], 'w');
    hold on
end
set(gca, 'YDir', 'Reverse'); 
set(gca, 'YTick', 1.1:1:6.1);
set(gca, 'YTickLabel', {'NW', 'CW', 'SW', 'SE', 'NE', 'All'});
xlabel('Percent');
axis([0 100 0.2 7]);
title('North Pacific Gyre Oscillation');
pbaspect([2 1 1]);
clear r l

clear ONI_Grid Other_Grid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Create maps of effort by phase
% Create lat and lon grids
% Horizontally stack column with latitudes
lat_grid = repmat(Lat,1,size(Lon,1));
% Transpose column with longitudes to a row and vertically stack
lon_grid = repmat(Lon',size(Lat,1),1);

% Filter effort to be non-confidential
% Permute so arrays are Lat x Lon x time
Vessels = permute(Vessels, [2 1 3]);
Effort_nonconfid = Effort;
Tot_vessels = sum(Vessels, 3, 'omitnan');
for r = 1:1:size(Tot_vessels,1)
    for c = 1:1:size(Tot_vessels,2)
        if Tot_vessels(r,c) < 3
            Effort_nonconfid(r,c,:) = 0;
        end
    end
end

% There's a function for this at the end of the script
% Effort_Map(effort_to_map, lat_grid_used, lon_grid_used, climate_mode_phase, panel_title)
% ONI
figure
subplot(2,2,1)
Effort_Map(Effort_nonconfid, lat_grid, lon_grid, PosONI, 'El Niño')

subplot(2,2,2)
Effort_Map(Effort_nonconfid, lat_grid, lon_grid, NegONI, 'La Niña')

subplot(2,2,3)
Effort_Map(Effort_nonconfid, lat_grid, lon_grid, NeuONI, 'Neutral')

subplot(2,2,4)
Effort_Map(Effort_nonconfid, lat_grid, lon_grid, NonONI, '<5 Consecutive Months')

% PDO
figure
subplot(1,2,1)
Effort_Map(Effort_nonconfid, lat_grid, lon_grid, PosPDO, 'Positive PDO')

subplot(1,2,2)
Effort_Map(Effort_nonconfid, lat_grid, lon_grid, NegPDO, 'Negative PDO')

% NPGO
figure
subplot(1,2,1)
Effort_Map(Effort_nonconfid, lat_grid, lon_grid, PosNPGO, 'Positive NPGO')

subplot(1,2,2)
Effort_Map(Effort_nonconfid, lat_grid, lon_grid, NegNPGO, 'Negative NPGO')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Functions
% Function to create effort + ONI plot (in case we want to make more)
function ONI_Effort_yyPlot(region_name, region_number, climate, fishery)
colororder({'#A5AAAf', '#737BE6'}) % y-axes colors, left to right
yyaxis left
plot(climate.ONI, 'Color', '#A5AAAf');
hold on
plot(climate.NeuONI, '-k', 'LineWidth', 2);
plot(climate.PosONI, '-r', 'LineWidth', 2);
plot(climate.NegONI, '-b', 'LineWidth', 2);
ylim([-3 3]);
ylabel('Oceanic Niño Index');
yyaxis right
plot(fishery(:,region_number), 'Color', '#737BE6');
ylim([0  6e6]);
ylabel('Hooks');
xlim([1 360]);
set(gca,'XTick',1:60:360);
set(gca,'XTickLabel',1995:5:2020);
title(sprintf('%s Region', region_name))
pbaspect([2 1 1]);
end


% Function to create effort + PDO plots
function PDO_Effort_yyPlot(region_name, region_number, climate, fishery)
colororder({'#A5AAAF', 'k'}) % y-axes colors, left to right
yyaxis left
stem(climate.PosPDO, '-r', 'Marker', 'none');
hold on
stem(climate.NegPDO, '-b', 'Marker', 'none');
ylim([-3 3]);
ylabel('Pacific Decadal Oscillation');
yyaxis right
plot(fishery(:,region_number));
ylim([0  6e6]);
ylabel('Hooks');
xlim([1 360]);
set(gca,'XTick',1:60:360);
set(gca,'XTickLabel',1995:5:2020);
title(sprintf('%s Region', region_name))
pbaspect([2 1 1]);
end


% Function to create effort + NPGO plots
function NPGO_Effort_yyPlot(region_name, region_number, climate, fishery)
colororder({'#A5AAAF', 'k'}) % y-axes colors, left to right
yyaxis left
stem(climate.PosNPGO, '-r', 'Marker', 'none');
hold on
stem(climate.NegNPGO, '-b', 'Marker', 'none');
ylim([-3.5 3.5]);
ylabel('North Pacific Gyre Oscillation');
yyaxis right
plot(fishery(:,region_number));
ylim([0  6e6]);
ylabel('Hooks');
xlim([1 360]);
set(gca,'XTick',1:60:360);
set(gca,'XTickLabel',1995:5:2020);
title(sprintf('%s Region', region_name))
pbaspect([2 1 1]);
end

% Function to create maps of effort by phase
function Effort_Map(effort_to_map, lat_grid_used, lon_grid_used, climate_mode_phase, panel_title)
axesm('mercator','MapLatLimit',[10 40],'MapLonLimit',[180 230], ...
    'MLineLocation', 10, 'PLineLocation', 10, ... % draw every 10 degrees         
    'Grid', 'on', 'MeridianLabel','on','ParallelLabel','on', ...
    'MLabelParallel', 10, 'MLabelLocation', 10, 'PLabelLocation', 10); % label every 10 degrees, below map 
Ig = geoshow(lat_grid_used,lon_grid_used,log10(sum(effort_to_map(:,:,~isnan(climate_mode_phase)),3)),'DisplayType','texturemap');
set(Ig,'AlphaData',double(~isinf(log10(sum(effort_to_map,3)))),'AlphaDataMapping','none','FaceAlpha','texturemap'); 
plotm([20 20], [180 230], 'k');
plotm([10 40], [210 210], 'k');
plotm([26 26], [180 210], 'k');
geoshow('landareas.shp','FaceColor',[0.5 0.5 0.5]);
clim([3 7])
cb = colorbar();
ylabel(cb, 'log_1_0(Effort)', 'Rotation', 270);
title(sprintf('%s'), panel_title)
set(gcf,'renderer','Painters')
tightmap
end













