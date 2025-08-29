% The purpose of this script is to look at any correlations between
% regional averages in oceanographic metrics and various modes of climate
% variability.

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
ONI = readtable('../ClimateIndices/ONI_withPhases.csv');
PDO = readtable('../ClimateIndices/PDO.csv');
NPGO = readtable('../ClimateIndices/NPGO.csv');
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Wrangle data
% Limit climate indices to our period of interest: 1995–2024
ONI_9524 = ONI(ONI.YR >= 1995 & ONI.YR <= 2024,:);
PDO_9524 = PDO(PDO.Year >= 1995 & PDO.Year <= 2024,:);
NPGO_9524 = NPGO(NPGO.YEAR >= 1995 & NPGO.YEAR <= 2024,:);
% Clean up
clear ONI PDO NPGO

% Permute so arrays are Lat x Lon (x Depth) x time
Catchability = permute(Catchability, [2 1 3]);
GODAS = permute(GODAS, [2 1 3 4]);
O2_2mlpl = permute(O2_2mlpl, [2 1 3]);

% Average various oceanographic metrics across regions
NW_lat = find(Lat >= 26);
CW_lat = find(Lat >= 20 & Lat < 26);
S_lat = find(Lat >= 10 & Lat < 20);
NE_lat = find(Lat >= 20);
W_lon = find(Lon >= 180 & Lon < 210);
E_lon = find (Lon >= 210 & Lon < 230);

Regional_Catchability(1:length(Month), 1:5) = NaN; %Cols are NW, CW, SW, SE, NE
Regional_O2depth(1:length(Month), 1:5) = NaN; %Cols are NW, CW, SW, SE, NE
Regional_Depth8deg(1:length(Month), 1:5) = NaN; %Cols are NW, CW, SW, SE, NE
Regional_Thickness814deg(1:length(Month), 1:5) = NaN; %Cols are NW, CW, SW, SE, NE
for m = 1:1:length(Month)
    % Catchability
    Regional_Catchability(m,1) = mean(Catchability(NW_lat, W_lon, m), "all", "omitnan");
    Regional_Catchability(m,2) = mean(Catchability(CW_lat, W_lon, m), "all", "omitnan");
    Regional_Catchability(m,3) = mean(Catchability(S_lat, W_lon, m), "all", "omitnan");
    Regional_Catchability(m,4) = mean(Catchability(S_lat, E_lon, m), "all", "omitnan");
    Regional_Catchability(m,5) = mean(Catchability(NE_lat, E_lon, m), "all", "omitnan");

    % Depth of 2 mlpl O2
    Regional_O2depth(m,1) = mean(O2_2mlpl(NW_lat, W_lon, m), "all", "omitnan");
    Regional_O2depth(m,2) = mean(O2_2mlpl(CW_lat, W_lon, m), "all", "omitnan");
    Regional_O2depth(m,3) = mean(O2_2mlpl(S_lat, W_lon, m), "all", "omitnan");
    Regional_O2depth(m,4) = mean(O2_2mlpl(S_lat, E_lon, m), "all", "omitnan");
    Regional_O2depth(m,5) = mean(O2_2mlpl(NE_lat, E_lon, m), "all", "omitnan");
    
    % For temperature, each month has r x c grid cells each with 31 depths
    % For each grid cell we're going to:
    % Find the max depth of 8-deg waters and
    % Determine the thickness of the 8-14-deg layer
    % Then, we'll average over space and fill the matrix
    Depth8(1:length(Lat), 1:length(Lon)) = NaN;
    Thick814(1:length(Lat), 1:length(Lon)) = NaN;

    for r = 1:1:length(Lat)
        for c = 1:1:length(Lon)
            Z8_idx = find(GODAS(r,c,:,m) >= 8);
            Z14_idx = find(GODAS(r,c,:,m) <= 14);
            
            Depth8(r,c) = GODAS_Depth(max(Z8_idx));
            Thick814(r,c) = GODAS_Depth(max(Z8_idx)) - GODAS_Depth(min(Z14_idx));
        end
    end
    clear r c *idx
    
    % Max depth of 8 deg C waters
    Regional_Depth8deg(m,1) = mean(Depth8(NW_lat, W_lon), "all", "omitnan");
    Regional_Depth8deg(m,2) = mean(Depth8(CW_lat, W_lon), "all", "omitnan");
    Regional_Depth8deg(m,3) = mean(Depth8(S_lat, W_lon), "all", "omitnan");
    Regional_Depth8deg(m,4) = mean(Depth8(S_lat, E_lon), "all", "omitnan");
    Regional_Depth8deg(m,5) = mean(Depth8(NE_lat, E_lon), "all", "omitnan");

    % Thickness of 8–14 deg C waters
    Regional_Thickness814deg(m,1) = mean(Thick814(NW_lat, W_lon), "all", "omitnan");
    Regional_Thickness814deg(m,2) = mean(Thick814(CW_lat, W_lon), "all", "omitnan");
    Regional_Thickness814deg(m,3) = mean(Thick814(S_lat, W_lon), "all", "omitnan");
    Regional_Thickness814deg(m,4) = mean(Thick814(S_lat, E_lon), "all", "omitnan");
    Regional_Thickness814deg(m,5) = mean(Thick814(NE_lat, E_lon), "all", "omitnan");
end
clear m NW_lat CW_lat S_lat NE_lat W_lon E_lon Depth8 Thick814

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Correlations
% Compare to climate variability index
% Empty arrays to fill (mode x region)
C_Spearman_r(1:3,1:5) = NaN;
C_Spearman_p(1:3,1:5) = NaN;
O_Spearman_r(1:3,1:5) = NaN;
O_Spearman_p(1:3,1:5) = NaN;
D_Spearman_r(1:3,1:5) = NaN;
D_Spearman_p(1:3,1:5) = NaN;
T_Spearman_r(1:3,1:5) = NaN;
T_Spearman_p(1:3,1:5) = NaN;

C_Pearson_r(1:3,1:5) = NaN;
C_Pearson_p(1:3,1:5) = NaN;
O_Pearson_r(1:3,1:5) = NaN;
O_Pearson_p(1:3,1:5) = NaN;
D_Pearson_r(1:3,1:5) = NaN;
D_Pearson_p(1:3,1:5) = NaN;
T_Pearson_r(1:3,1:5) = NaN;
T_Pearson_p(1:3,1:5) = NaN;

for r = 1:1:5
    % ONI
    [C_Spearman_r(1,r), C_Spearman_p(1,r)] = corr(Regional_Catchability(:,r), ONI_9524.ONI,...
        'type', 'Spearman', 'rows', 'complete');
    [O_Spearman_r(1,r), O_Spearman_p(1,r)] = corr(Regional_O2depth(:,r), ONI_9524.ONI,...
        'type', 'Spearman', 'rows', 'complete');
    [D_Spearman_r(1,r), D_Spearman_p(1,r)] = corr(Regional_Depth8deg(:,r), ONI_9524.ONI,...
        'type', 'Spearman', 'rows', 'complete');
    [T_Spearman_r(1,r), T_Spearman_p(1,r)] = corr(Regional_Thickness814deg(:,r), ONI_9524.ONI,...
        'type', 'Spearman', 'rows', 'complete');

    [Cr1, Cp1] = corrcoef(Regional_Catchability(:,r), ONI_9524.ONI, 'Rows', 'pairwise');
    [Or1, Op1] = corrcoef(Regional_O2depth(:,r), ONI_9524.ONI, 'Rows', 'pairwise');
    [Dr1, Dp1] = corrcoef(Regional_Depth8deg(:,r), ONI_9524.ONI, 'Rows', 'pairwise');
    [Tr1, Tp1] = corrcoef(Regional_Thickness814deg(:,r), ONI_9524.ONI, 'Rows', 'pairwise');

    C_Pearson_r(1,r) = Cr1(1,2);
    C_Pearson_p(1,r) = Cp1(1,2);
    O_Pearson_r(1,r) = Or1(1,2);
    O_Pearson_p(1,r) = Op1(1,2);
    D_Pearson_r(1,r) = Dr1(1,2);
    D_Pearson_p(1,r) = Dp1(1,2);
    T_Pearson_r(1,r) = Tr1(1,2);
    T_Pearson_p(1,r) = Tp1(1,2);

    % PDO  
    [C_Spearman_r(2,r), C_Spearman_p(2,r)] = corr(Regional_Catchability(:,r), PDO_9524.PDO,...
        'type', 'Spearman', 'rows', 'complete');
    [O_Spearman_r(2,r), O_Spearman_p(2,r)] = corr(Regional_O2depth(:,r), PDO_9524.PDO,...
        'type', 'Spearman', 'rows', 'complete');
    [D_Spearman_r(2,r), D_Spearman_p(2,r)] = corr(Regional_Depth8deg(:,r), PDO_9524.PDO,...
        'type', 'Spearman', 'rows', 'complete');
    [T_Spearman_r(2,r), T_Spearman_p(2,r)] = corr(Regional_Thickness814deg(:,r), PDO_9524.PDO,...
        'type', 'Spearman', 'rows', 'complete');

    [Cr2, Cp2] = corrcoef(Regional_Catchability(:,r), PDO_9524.PDO, 'Rows', 'pairwise');
    [Or2, Op2] = corrcoef(Regional_O2depth(:,r), PDO_9524.PDO, 'Rows', 'pairwise');
    [Dr2, Dp2] = corrcoef(Regional_Depth8deg(:,r), PDO_9524.PDO, 'Rows', 'pairwise');
    [Tr2, Tp2] = corrcoef(Regional_Thickness814deg(:,r), PDO_9524.PDO, 'Rows', 'pairwise');

    C_Pearson_r(2,r) = Cr2(1,2);
    C_Pearson_p(2,r) = Cp2(1,2);
    O_Pearson_r(2,r) = Or2(1,2);
    O_Pearson_p(2,r) = Op2(1,2);
    D_Pearson_r(2,r) = Dr2(1,2);
    D_Pearson_p(2,r) = Dp2(1,2);
    T_Pearson_r(2,r) = Tr2(1,2);
    T_Pearson_p(2,r) = Tp2(1,2);

    % NPGO 
    [C_Spearman_r(3,r), C_Spearman_p(3,r)] = corr(Regional_Catchability(:,r),NPGO_9524.NPGO,...
        'type', 'Spearman', 'rows', 'complete');
    [O_Spearman_r(3,r), O_Spearman_p(3,r)] = corr(Regional_O2depth(:,r), NPGO_9524.NPGO,...
        'type', 'Spearman', 'rows', 'complete');
    [D_Spearman_r(3,r), D_Spearman_p(3,r)] = corr(Regional_Depth8deg(:,r), NPGO_9524.NPGO,...
        'type', 'Spearman', 'rows', 'complete');
    [T_Spearman_r(3,r), T_Spearman_p(3,r)] = corr(Regional_Thickness814deg(:,r), NPGO_9524.NPGO,...
        'type', 'Spearman', 'rows', 'complete');

    [Cr3, Cp3] = corrcoef(Regional_Catchability(:,r), NPGO_9524.NPGO, 'Rows', 'pairwise');
    [Or3, Op3] = corrcoef(Regional_O2depth(:,r), NPGO_9524.NPGO, 'Rows', 'pairwise');
    [Dr3, Dp3] = corrcoef(Regional_Depth8deg(:,r), NPGO_9524.NPGO, 'Rows', 'pairwise');
    [Tr3, Tp3] = corrcoef(Regional_Thickness814deg(:,r), NPGO_9524.NPGO, 'Rows', 'pairwise');

    C_Pearson_r(3,r) = Cr3(1,2);
    C_Pearson_p(3,r) = Cp3(1,2);
    O_Pearson_r(3,r) = Or3(1,2);
    O_Pearson_p(3,r) = Op3(1,2);
    D_Pearson_r(3,r) = Dr3(1,2);
    D_Pearson_p(3,r) = Dp3(1,2);
    T_Pearson_r(3,r) = Tr3(1,2);
    T_Pearson_p(3,r) = Tp3(1,2);
end
clear r

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plots 
% Plots for significant correlations
% There are functions for this at the end
% Catchability_Plot(ocean_data, mode_name, mode_data, regions)
figure
Catchability_Plot(Regional_Catchability, "ONI", ONI_9524.ONI, 1:1:4)
% This proved busy and not terribly informative, so I stopped here

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Functions
function Catchability_Plot(ocean_data, mode_name, mode_data, regions)
RGB = orderedcolors("gem"); % colors for region (Matlab defaults)
colororder({'#A5AAAf', 'k'}); % y-axes colors, left (grey) to right (black)
yyaxis left
for r = 1:1:length(regions)
    plot(1:360, ocean_data(:,r), 'Color', RGB(r,:), 'LineStyle', '-');
    hold on
end
ylim([0 0.6]);
ylabel('Catchability');
yyaxis right
plot(1:360, mode_data, 'k');
ylim([-3.5 3.5]);
ylabel('Index Value');
xlim([1 360]);
set(gca,'XTick',1:60:360);
set(gca,'XTickLabel',1995:5:2020);
title(sprintf('%s Region', mode_name));
pbaspect([2 1 1]);
end

















