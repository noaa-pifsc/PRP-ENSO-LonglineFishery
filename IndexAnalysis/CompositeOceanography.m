% The purpose of this script is to look at composite oceanographic states
% for each phase of each mode of climate variability.
% We'll also look at spatiatemploral correlations because much of the prep
% work is the same.

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

% To make composites, it is helpful to separate out positive and negative
% phases of each index (which is a little tedious)
PosONI = ONI_9524.ONI;
PosONI(isnan(str2double(ONI_9524.ElNino))) = NaN;
NegONI = ONI_9524.ONI;
NegONI(isnan(str2double(ONI_9524.LaNina))) = NaN;
NeuONI = ONI_9524.ONI;
NeuONI(isnan(str2double(ONI_9524.Neutral))) = NaN;
NonONI = ONI_9524.ONI;
NonONI(~isnan(str2double(ONI_9524.ElNino)) | ~isnan(str2double(ONI_9524.LaNina)) | ~isnan(str2double(ONI_9524.Neutral))) = NaN;

PosPDO = PDO_9524.PDO;
PosPDO(PosPDO < 0) = NaN;
NegPDO = PDO_9524.PDO;
NegPDO(NegPDO > 0) = NaN;
NeuPDO = PDO_9524.PDO;
NeuPDO(NeuPDO ~= 0) = NaN;

PosNPGO = NPGO_9524.NPGO;
PosNPGO(PosNPGO < 0) = NaN;
NegNPGO = NPGO_9524.NPGO;
NegNPGO(NegNPGO > 0) = NaN;
NeuNPGO = NPGO_9524.NPGO;
NeuNPGO(NeuNPGO ~= 0) = NaN;

% Permute so arrays are Lat x Lon (x Depth) x time
Catchability = permute(Catchability, [2 1 3]);
GODAS = permute(GODAS, [2 1 3 4]);
O2_2mlpl = permute(O2_2mlpl, [2 1 3]);

% Create lat and lon grids (for plotting)
% Horizontally stack column with latitudes
lat_grid = repmat(Lat,1,size(Lon,1));
% Transpose column with longitudes to a row and vertically stack
lon_grid = repmat(Lon',size(Lat,1),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Composite states
% One way to approach this is to look at the difference from the mean state
Catch_climo = mean(Catchability, 3, "omitnan");
% GODAS_climo = mean(GODAS, 4, "omitnan"); % I'm not sure this is the best approach here
O2_climo = mean(O2_2mlpl, 3, "omitnan");

PosONI_Catch_diff = mean(Catchability(:,:,~isnan(PosONI)), 3, "omitnan") - Catch_climo;
NegONI_Catch_diff = mean(Catchability(:,:,~isnan(NegONI)), 3, "omitnan") - Catch_climo;
NeuONI_Catch_diff = mean(Catchability(:,:,~isnan(NeuONI)), 3, "omitnan") - Catch_climo;
NonONI_Catch_diff = mean(Catchability(:,:,~isnan(NonONI)), 3, "omitnan") - Catch_climo;

PosPDO_Catch_diff = mean(Catchability(:,:,~isnan(PosPDO)), 3, "omitnan") - Catch_climo;
NegPDO_Catch_diff = mean(Catchability(:,:,~isnan(NegPDO)), 3, "omitnan") - Catch_climo;

PosNPGO_Catch_diff = mean(Catchability(:,:,~isnan(PosNPGO)), 3, "omitnan") - Catch_climo;
NegNPGO_Catch_diff = mean(Catchability(:,:,~isnan(NegNPGO)), 3, "omitnan") - Catch_climo;

PosONI_O2_diff = mean(O2_2mlpl(:,:,~isnan(PosONI)), 3, "omitnan") - O2_climo;
NegONI_O2_diff = mean(O2_2mlpl(:,:,~isnan(NegONI)), 3, "omitnan") - O2_climo;
NeuONI_O2_diff = mean(O2_2mlpl(:,:,~isnan(NeuONI)), 3, "omitnan") - O2_climo;
NonONI_O2_diff = mean(O2_2mlpl(:,:,~isnan(NonONI)), 3, "omitnan") - O2_climo;

PosPDO_O2_diff = mean(O2_2mlpl(:,:,~isnan(PosPDO)), 3, "omitnan") - O2_climo;
NegPDO_O2_diff = mean(O2_2mlpl(:,:,~isnan(NegPDO)), 3, "omitnan") - O2_climo;

PosNPGO_O2_diff = mean(O2_2mlpl(:,:,~isnan(PosNPGO)), 3, "omitnan") - O2_climo;
NegNPGO_O2_diff = mean(O2_2mlpl(:,:,~isnan(NegNPGO)), 3, "omitnan") - O2_climo;

% One approach to temperature:
% Max depth of 8-deg C and thickness of 8-14-deg C layer
% Empty matrices to fill
GODAS_8degDepth_climo(1:length(Lat), 1:length(Lon)) = NaN;
GODAS_8to14degThickness_climo(1:length(Lat), 1:length(Lon)) = NaN;

PosONI_GODAS_8degDepth_diff(1:length(Lat), 1:length(Lon)) = NaN;
NegONI_GODAS_8degDepth_diff(1:length(Lat), 1:length(Lon)) = NaN;
NeuONI_GODAS_8degDepth_diff(1:length(Lat), 1:length(Lon)) = NaN;
NonONI_GODAS_8degDepth_diff(1:length(Lat), 1:length(Lon)) = NaN;

PosPDO_GODAS_8degDepth_diff(1:length(Lat), 1:length(Lon)) = NaN;
NegPDO_GODAS_8degDepth_diff(1:length(Lat), 1:length(Lon)) = NaN;

PosNPGO_GODAS_8degDepth_diff(1:length(Lat), 1:length(Lon)) = NaN;
NegNPGO_GODAS_8degDepth_diff(1:length(Lat), 1:length(Lon)) = NaN;

PosONI_GODAS_8to14degThickness_diff(1:length(Lat), 1:length(Lon)) = NaN;
NegONI_GODAS_8to14degThickness_diff(1:length(Lat), 1:length(Lon)) = NaN;
NeuONI_GODAS_8to14degThickness_diff(1:length(Lat), 1:length(Lon)) = NaN;
NonONI_GODAS_8to14degThickness_diff(1:length(Lat), 1:length(Lon)) = NaN;

PosPDO_GODAS_8to14degThickness_diff(1:length(Lat), 1:length(Lon)) = NaN;
NegPDO_GODAS_8to14degThickness_diff(1:length(Lat), 1:length(Lon)) = NaN;

PosNPGO_GODAS_8to14degThickness_diff(1:length(Lat), 1:length(Lon)) = NaN;
NegNPGO_GODAS_8to14degThickness_diff(1:length(Lat), 1:length(Lon)) = NaN;

for r = 1:1:length(Lat)
    for c = 1:1:length(Lon)
        % Each grid cell has 31 depths at each of 360 months
        % At each time step, we're going to:
        % Find the max depth of 8-deg waters and
        % Determine the thickness of the 8-14-deg layer
        % Then, we'll average over time and fill the matrix
        Depth8(1:360) = NaN;
        Thick814(1:360) = NaN;
        for m = 1:1:360
            Z8_idx = find(GODAS(r,c,:,m) >= 8);
            Z14_idx = find(GODAS(r,c,:,m) <= 14);
            
            Depth8(m) = GODAS_Depth(max(Z8_idx));
            Thick814(m) = GODAS_Depth(max(Z8_idx)) - GODAS_Depth(min(Z14_idx));
            
        end
        clear m *idx

        GODAS_8degDepth_climo(r,c) = mean(Depth8, "omitnan");
        GODAS_8to14degThickness_climo(r,c) = mean(Thick814, "omitnan");

        % We'll go ahead and do each of the phases of each of the modes of
        % variability, too.
        PosONI_GODAS_8degDepth_diff(r,c) = mean(Depth8(~isnan(PosONI)), "omitnan") - mean(Depth8, "omitnan");
        NegONI_GODAS_8degDepth_diff(r,c) = mean(Depth8(~isnan(NegONI)), "omitnan") - mean(Depth8, "omitnan");
        NeuONI_GODAS_8degDepth_diff(r,c) = mean(Depth8(~isnan(NeuONI)), "omitnan") - mean(Depth8, "omitnan");
        NonONI_GODAS_8degDepth_diff(r,c) = mean(Depth8(~isnan(NonONI)), "omitnan") - mean(Depth8, "omitnan");

        PosPDO_GODAS_8degDepth_diff(r,c) = mean(Depth8(~isnan(PosPDO)), "omitnan") - mean(Depth8, "omitnan");
        NegPDO_GODAS_8degDepth_diff(r,c) = mean(Depth8(~isnan(NegPDO)), "omitnan") - mean(Depth8, "omitnan");

        PosNPGO_GODAS_8degDepth_diff(r,c) = mean(Depth8(~isnan(PosNPGO)), "omitnan") - mean(Depth8, "omitnan");
        NegNPGO_GODAS_8degDepth_diff(r,c) = mean(Depth8(~isnan(NegNPGO)), "omitnan") - mean(Depth8, "omitnan");

        PosONI_GODAS_8to14degThickness_diff(r,c) = mean(Thick814(~isnan(PosONI)), "omitnan") - mean(Thick814, "omitnan");
        NegONI_GODAS_8to14degThickness_diff(r,c) = mean(Thick814(~isnan(NegONI)), "omitnan") - mean(Thick814, "omitnan");
        NeuONI_GODAS_8to14degThickness_diff(r,c) = mean(Thick814(~isnan(NeuONI)), "omitnan") - mean(Thick814, "omitnan");
        NonONI_GODAS_8to14degThickness_diff(r,c) = mean(Thick814(~isnan(NonONI)), "omitnan") - mean(Thick814, "omitnan");

        PosPDO_GODAS_8to14degThickness_diff(r,c) = mean(Thick814(~isnan(PosPDO)), "omitnan") - mean(Thick814, "omitnan");
        NegPDO_GODAS_8to14degThickness_diff(r,c) = mean(Thick814(~isnan(NegPDO)), "omitnan") - mean(Thick814, "omitnan");
        
        PosNPGO_GODAS_8to14degThickness_diff(r,c) = mean(Thick814(~isnan(PosNPGO)), "omitnan") - mean(Thick814, "omitnan");
        NegNPGO_GODAS_8to14degThickness_diff(r,c) = mean(Thick814(~isnan(NegNPGO)), "omitnan") - mean(Thick814, "omitnan");

    end
end
clear r c Depth8 Thick814

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Correlations
% Empty array to fill (lat x lon x mode: ENSO, PDO, NPGO)
Catch_Spearman_r(1:length(Lat), 1:length(Lon), 1:3) = NaN;
Catch_Spearman_p(1:length(Lat), 1:length(Lon), 1:3) = NaN;
Catch_Pearson_r(1:length(Lat), 1:length(Lon), 1:3) = NaN;
Catch_Pearson_p(1:length(Lat), 1:length(Lon), 1:3) = NaN;

O2_Spearman_r(1:length(Lat), 1:length(Lon), 1:3) = NaN;
O2_Spearman_p(1:length(Lat), 1:length(Lon), 1:3) = NaN;
O2_Pearson_r(1:length(Lat), 1:length(Lon), 1:3) = NaN;
O2_Pearson_p(1:length(Lat), 1:length(Lon), 1:3) = NaN;

D8_Spearman_r(1:length(Lat), 1:length(Lon), 1:3) = NaN;
D8_Spearman_p(1:length(Lat), 1:length(Lon), 1:3) = NaN;
D8_Pearson_r(1:length(Lat), 1:length(Lon), 1:3) = NaN;
D8_Pearson_p(1:length(Lat), 1:length(Lon), 1:3) = NaN;

T814_Spearman_r(1:length(Lat), 1:length(Lon), 1:3) = NaN;
T814_Spearman_p(1:length(Lat), 1:length(Lon), 1:3) = NaN;
T814_Pearson_r(1:length(Lat), 1:length(Lon), 1:3) = NaN;
T814_Pearson_p(1:length(Lat), 1:length(Lon), 1:3) = NaN;

for r = 1:1:length(Lat)
    for c = 1:1:length(Lon)
        Catch = squeeze(Catchability(r, c, :));
        O2 = squeeze(O2_2mlpl(r, c, :));

        Depth8(1:360) = NaN;
        Thick814(1:360) = NaN;
        for m = 1:1:360
            Z8_idx = find(GODAS(r,c,:,m) >= 8);
            Z14_idx = find(GODAS(r,c,:,m) <= 14);
            
            Depth8(m) = GODAS_Depth(max(Z8_idx));
            Thick814(m) = GODAS_Depth(max(Z8_idx)) - GODAS_Depth(min(Z14_idx));
            
        end
        clear m *idx 

        % Catchability with index value
        [Catch_Spearman_r(r, c, 1), Catch_Spearman_p(r, c, 1)] = ...
            corr(Catch, ONI_9524.ONI, 'type', 'Spearman', 'rows', 'complete');
        [Catch_Spearman_r(r, c, 2), Catch_Spearman_p(r, c, 2)] = ...
            corr(Catch, PDO_9524.PDO, 'type', 'Spearman', 'rows', 'complete');
        [Catch_Spearman_r(r, c, 3), Catch_Spearman_p(r, c, 3)] = ...
            corr(Catch, NPGO_9524.NPGO, 'type', 'Spearman', 'rows', 'complete');

        [CPr1, CPp1] = corrcoef(Catch, ONI_9524.ONI, 'Rows', 'pairwise');
        [CPr2, CPp2] = corrcoef(Catch, PDO_9524.PDO, 'Rows', 'pairwise');
        [CPr3, CPp3] = corrcoef(Catch, NPGO_9524.NPGO, 'Rows', 'pairwise');

        Catch_Pearson_r(r, c, 1) = CPr1(1,2);
        Catch_Pearson_p(r, c, 1) = CPp1(1,2);
        Catch_Pearson_r(r, c, 2) = CPr2(1,2);
        Catch_Pearson_p(r, c, 2) = CPp2(1,2);
        Catch_Pearson_r(r, c, 3) = CPr3(1,2);
        Catch_Pearson_p(r, c, 3) = CPp3(1,2);

        % Depth of 2mlpl O2 with index value
        [O2_Spearman_r(r, c, 1), O2_Spearman_p(r, c, 1)] = ...
            corr(O2, ONI_9524.ONI, 'type', 'Spearman', 'rows', 'complete');
        [O2_Spearman_r(r, c, 2), O2_Spearman_p(r, c, 2)] = ...
            corr(O2, PDO_9524.PDO, 'type', 'Spearman', 'rows', 'complete');
        [O2_Spearman_r(r, c, 3), O2_Spearman_p(r, c, 3)] = ...
            corr(O2, NPGO_9524.NPGO, 'type', 'Spearman', 'rows', 'complete');

        [OPr1, OPp1] = corrcoef(O2, ONI_9524.ONI, 'Rows', 'pairwise');
        [OPr2, OPp2] = corrcoef(O2, PDO_9524.PDO, 'Rows', 'pairwise');
        [OPr3, OPp3] = corrcoef(O2, NPGO_9524.NPGO, 'Rows', 'pairwise');

        O2_Pearson_r(r, c, 1) = OPr1(1,2);
        O2_Pearson_p(r, c, 1) = OPp1(1,2);
        O2_Pearson_r(r, c, 2) = OPr2(1,2);
        O2_Pearson_p(r, c, 2) = OPp2(1,2);
        O2_Pearson_r(r, c, 3) = OPr3(1,2);
        O2_Pearson_p(r, c, 3) = OPp3(1,2);

        % Depth of 8 deg C water with index value
        [D8_Spearman_r(r, c, 1), D8_Spearman_p(r, c, 1)] = ...
            corr(Depth8', ONI_9524.ONI, 'type', 'Spearman', 'rows', 'complete');
        [D8_Spearman_r(r, c, 2), D8_Spearman_p(r, c, 2)] = ...
            corr(Depth8', PDO_9524.PDO, 'type', 'Spearman', 'rows', 'complete');
        [D8_Spearman_r(r, c, 3), D8_Spearman_p(r, c, 3)] = ...
            corr(Depth8', NPGO_9524.NPGO, 'type', 'Spearman', 'rows', 'complete');

        [DPr1, DPp1] = corrcoef(Depth8, ONI_9524.ONI, 'Rows', 'pairwise');
        [DPr2, DPp2] = corrcoef(Depth8, PDO_9524.PDO, 'Rows', 'pairwise');
        [DPr3, DPp3] = corrcoef(Depth8, NPGO_9524.NPGO, 'Rows', 'pairwise');

        D8_Pearson_r(r, c, 1) = DPr1(1,2);
        D8_Pearson_p(r, c, 1) = DPp1(1,2);
        D8_Pearson_r(r, c, 2) = DPr2(1,2);
        D8_Pearson_p(r, c, 2) = DPp2(1,2);
        D8_Pearson_r(r, c, 3) = DPr3(1,2);
        D8_Pearson_p(r, c, 3) = DPp3(1,2);

        % Thickness of 8-14 deg C water with index value
        [T814_Spearman_r(r, c, 1), T814_Spearman_p(r, c, 1)] = ...
            corr(Thick814', ONI_9524.ONI, 'type', 'Spearman', 'rows', 'complete');
        [T814_Spearman_r(r, c, 2), T814_Spearman_p(r, c, 2)] = ...
            corr(Thick814', PDO_9524.PDO, 'type', 'Spearman', 'rows', 'complete');
        [T814_Spearman_r(r, c, 3), T814_Spearman_p(r, c, 3)] = ...
            corr(Thick814', NPGO_9524.NPGO, 'type', 'Spearman', 'rows', 'complete');

        [TPr1, TPp1] = corrcoef(Thick814, ONI_9524.ONI, 'Rows', 'pairwise');
        [TPr2, TPp2] = corrcoef(Thick814, PDO_9524.PDO, 'Rows', 'pairwise');
        [TPr3, TPp3] = corrcoef(Thick814, NPGO_9524.NPGO, 'Rows', 'pairwise');

        T814_Pearson_r(r, c, 1) = TPr1(1,2);
        T814_Pearson_p(r, c, 1) = TPp1(1,2);
        T814_Pearson_r(r, c, 2) = TPr2(1,2);
        T814_Pearson_p(r, c, 2) = TPp2(1,2);
        T814_Pearson_r(r, c, 3) = TPr3(1,2);
        T814_Pearson_p(r, c, 3) = TPp3(1,2);

    end
end
clear r c CP* DP* OP* Depth8 Thick814 Catch O2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plots
% There are functions for this at the end of the script

% Catchability
% figure
% Composite_Map(lat_grid, lon_grid, Catch_climo, 0, 1, "Catchability Climatology")
% 
% figure
% subplot(2,2,1)
% Difference_Map(lat_grid, lon_grid, PosONI_Catch_diff, -0.06, 0.06, "Catchability: El Niño")
% subplot(2,2,2)
% Difference_Map(lat_grid, lon_grid, NegONI_Catch_diff, -0.06, 0.06, "Catchability: La Niña")
% subplot(2,2,3)
% Difference_Map(lat_grid, lon_grid, NeuONI_Catch_diff, -0.06, 0.06, "Catchability: ENSO-Neutral")
% subplot(2,2,4)
% Difference_Map(lat_grid, lon_grid, NonONI_Catch_diff, -0.06, 0.06, "Catchability: ENSO-No Phase")
% 
% figure
% subplot(1,2,1)
% Difference_Map(lat_grid, lon_grid, PosPDO_Catch_diff, -0.06, 0.06, "Catchability: Positive PDO")
% subplot(1,2,2)
% Difference_Map(lat_grid, lon_grid, NegPDO_Catch_diff, -0.06, 0.06, "Catchability: Negative PDO")
% 
% figure
% subplot(1,2,1)
% Difference_Map(lat_grid, lon_grid, PosNPGO_Catch_diff, -0.06, 0.06, "Catchability: Positive NPGO")
% subplot(1,2,2)
% Difference_Map(lat_grid, lon_grid, NegNPGO_Catch_diff, -0.06, 0.06, "Catchability: Negative NPGO")
%
% figure
% subplot(1,2,1)
% Corr_Map(lat_grid, lon_grid, Catch_Spearman_r(:,:,1), Catch_Spearman_p(:,:,1), ...
%     0.05, -1, 1, "Spearman Correlation: ONI and Catchability")
% subplot(1,2,2)
% Corr_Map(lat_grid, lon_grid, Catch_Pearson_r(:,:,1), Catch_Pearson_p(:,:,1), ...
%     0.05, -1, 1, "Pearson Correlation: ONI and Catchability")
% 
% figure
% subplot(1,2,1)
% Corr_Map(lat_grid, lon_grid, Catch_Spearman_r(:,:,2), Catch_Spearman_p(:,:,2), ...
%     0.05, -1, 1, "Spearman Correlation: PDO and Catchability")
% subplot(1,2,2)
% Corr_Map(lat_grid, lon_grid, Catch_Pearson_r(:,:,2), Catch_Pearson_p(:,:,2), ...
%     0.05, -1, 1, "Pearson Correlation: PDO and Catchability")
% 
% figure
% subplot(1,2,1)
% Corr_Map(lat_grid, lon_grid, Catch_Spearman_r(:,:,3), Catch_Spearman_p(:,:,3), ...
%     0.05, -1, 1, "Spearman Correlation: NPGO and Catchability")
% subplot(1,2,2)
% Corr_Map(lat_grid, lon_grid, Catch_Pearson_r(:,:,3), Catch_Pearson_p(:,:,3), ...
%     0.05, -1, 1, "Pearson Correlation: NPGO and Catchability")

% Oxygen
% figure
% Composite_Map(lat_grid, lon_grid, O2_climo, 0, 800, "Climatology: depth of 2 mlpl O_2")
% 
% figure
% subplot(2,2,1)
% Difference_Map(lat_grid, lon_grid, PosONI_O2_diff, -25, 25, "Depth of 2 mlpl O_2: El Niño")
% subplot(2,2,2)
% Difference_Map(lat_grid, lon_grid, NegONI_O2_diff, -25, 25, "Depth of 2 mlpl O_2: La Niña")
% subplot(2,2,3)
% Difference_Map(lat_grid, lon_grid, NeuONI_O2_diff, -25, 25, "Depth of 2 mlpl O_2: ENSO-Neutral")
% subplot(2,2,4)
% Difference_Map(lat_grid, lon_grid, NonONI_O2_diff, -25, 25, "Depth of 2 mlpl O_2: ENSO-No Phase")
% 
% figure
% subplot(1,2,1)
% Difference_Map(lat_grid, lon_grid, PosPDO_O2_diff, -25, 25, "Depth of 2 mlpl O_2: Positive PDO")
% subplot(1,2,2)
% Difference_Map(lat_grid, lon_grid, NegPDO_O2_diff, -25, 25, "Depth of 2 mlpl O_2: Negative PDO")
% 
% figure
% subplot(1,2,1)
% Difference_Map(lat_grid, lon_grid, PosNPGO_O2_diff, -25, 25, "Depth of 2 mlpl O_2: Positive NPGO")
% subplot(1,2,2)
% Difference_Map(lat_grid, lon_grid, NegNPGO_O2_diff, -25, 25, "Depth of 2 mlpl O_2: Negative NPGO")
% 
% figure
% subplot(1,2,1)
% Corr_Map(lat_grid, lon_grid, O2_Spearman_r(:,:,1), O2_Spearman_p(:,:,1), ...
%     0.05, -1, 1, "Spearman Correlation: ONI and Depth of 2 mlpl O_2")
% subplot(1,2,2)
% Corr_Map(lat_grid, lon_grid, O2_Pearson_r(:,:,1), O2_Pearson_p(:,:,1), ...
%     0.05, -1, 1, "Pearson Correlation: ONI and Depth of 2 mlpl O_2")
% 
% figure
% subplot(1,2,1)
% Corr_Map(lat_grid, lon_grid, O2_Spearman_r(:,:,2), O2_Spearman_p(:,:,2), ...
%     0.05, -1, 1, "Spearman Correlation: PDO and Depth of 2 mlpl O_2")
% subplot(1,2,2)
% Corr_Map(lat_grid, lon_grid, O2_Pearson_r(:,:,2), O2_Pearson_p(:,:,2), ...
%     0.05, -1, 1, "Pearson Correlation: PDO and Depth of 2 mlpl O_2")
% 
% figure
% subplot(1,2,1)
% Corr_Map(lat_grid, lon_grid, O2_Spearman_r(:,:,3), O2_Spearman_p(:,:,3), ...
%     0.05, -1, 1, "Spearman Correlation: NPGO and Depth of 2 mlpl O_2")
% subplot(1,2,2)
% Corr_Map(lat_grid, lon_grid, O2_Pearson_r(:,:,3), O2_Pearson_p(:,:,3), ...
%     0.05, -1, 1, "Pearson Correlation: NPGO and Depth of 2 mlpl O_2")

% 8-14 deg C waters
figure
Thickness_Map(lat_grid, lon_grid, GODAS_8degDepth_climo, GODAS_8to14degThickness_climo, ...
    0, 700, "Climatology: Max depth (shaded) and thickness (contoured) of 8–14 deg C")

figure
subplot(2,2,1)
Thickness_Diff_Map(lat_grid, lon_grid, PosONI_GODAS_8degDepth_diff, PosONI_GODAS_8to14degThickness_diff, ...
    -25, 25, "El Niño: Max depth (shaded) and thickness (contoured) of 8–14 deg C")
subplot(2,2,2)
Thickness_Diff_Map(lat_grid, lon_grid, NegONI_GODAS_8degDepth_diff, NegONI_GODAS_8to14degThickness_diff, ...
    -25, 25, "La Niña: Max depth (shaded) and thickness (contoured) of 8–14 deg C")
subplot(2,2,3)
Thickness_Diff_Map(lat_grid, lon_grid, NeuONI_GODAS_8degDepth_diff, NeuONI_GODAS_8to14degThickness_diff, ...
    -25, 25, "ENSO Neutral: Max depth (shaded) and thickness (contoured) of 8–14 deg C")
subplot(2,2,4)
Thickness_Diff_Map(lat_grid, lon_grid, NonONI_GODAS_8degDepth_diff, NonONI_GODAS_8to14degThickness_diff, ...
    -25, 25, "ENSO No Phase: Max depth (shaded) and thickness (contoured) of 8–14 deg C")

figure
subplot(1,2,1)
Thickness_Diff_Map(lat_grid, lon_grid, PosPDO_GODAS_8degDepth_diff, PosPDO_GODAS_8to14degThickness_diff, ...
    -25, 25, "Positive PDO: Max depth (shaded) and thickness (contoured) of 8–14 deg C")
subplot(1,2,2)
Thickness_Diff_Map(lat_grid, lon_grid, NegPDO_GODAS_8degDepth_diff, NegPDO_GODAS_8to14degThickness_diff, ...
    -25, 25, "Negative PDO: Max depth (shaded) and thickness (contoured) of 8–14 deg C")

figure
subplot(1,2,1)
Thickness_Diff_Map(lat_grid, lon_grid, PosNPGO_GODAS_8degDepth_diff, PosNPGO_GODAS_8to14degThickness_diff, ...
    -25, 25, "Positive NPGO: Max depth (shaded) and thickness (contoured) of 8–14 deg C")
subplot(1,2,2)
Thickness_Diff_Map(lat_grid, lon_grid, NegNPGO_GODAS_8degDepth_diff, NegNPGO_GODAS_8to14degThickness_diff, ...
    -25, 25, "Negative NPGO: Max depth (shaded) and thickness (contoured) of 8–14 deg C")

figure
subplot(2,2,1)
Corr_Map(lat_grid, lon_grid, D8_Spearman_r(:,:,1), D8_Spearman_p(:,:,1), ...
    0.05, -1, 1, "Spearman Correlation: ONI and Max Depth of 8 C Temp")
subplot(2,2,2)
Corr_Map(lat_grid, lon_grid, D8_Pearson_r(:,:,1), D8_Pearson_p(:,:,1), ...
    0.05, -1, 1, "Pearson Correlation: ONI and Max Depth of 8 C Temp")
subplot(2,2,3)
Corr_Map(lat_grid, lon_grid, T814_Spearman_r(:,:,1), T814_Spearman_p(:,:,1), ...
    0.05, -1, 1, "Spearman Correlation: ONI and Thickness of 8-14 C Temps")
subplot(2,2,4)
Corr_Map(lat_grid, lon_grid, T814_Pearson_r(:,:,1), T814_Pearson_p(:,:,1), ...
    0.05, -1, 1, "Pearson Correlation: ONI and Thickness of 8-14 C Temps")

figure
subplot(2,2,1)
Corr_Map(lat_grid, lon_grid, D8_Spearman_r(:,:,2), D8_Spearman_p(:,:,2), ...
    0.05, -1, 1, "Spearman Correlation: PDO and Max Depth of 8 C Temp")
subplot(2,2,2)
Corr_Map(lat_grid, lon_grid, D8_Pearson_r(:,:,2), D8_Pearson_p(:,:,2), ...
    0.05, -1, 1, "Pearson Correlation: PDO and Max Depth of 8 C Temp")
subplot(2,2,3)
Corr_Map(lat_grid, lon_grid, T814_Spearman_r(:,:,2), T814_Spearman_p(:,:,2), ...
    0.05, -1, 1, "Spearman Correlation: PDO and Thickness of 8-14 C Temps")
subplot(2,2,4)
Corr_Map(lat_grid, lon_grid, T814_Pearson_r(:,:,2), T814_Pearson_p(:,:,2), ...
    0.05, -1, 1, "Pearson Correlation: PDO and Thickness of 8-14 C Temps")

figure
subplot(2,2,1)
Corr_Map(lat_grid, lon_grid, D8_Spearman_r(:,:,3), D8_Spearman_p(:,:,3), ...
    0.05, -1, 1, "Spearman Correlation: NPGO and Max Depth of 8 C Temp")
subplot(2,2,2)
Corr_Map(lat_grid, lon_grid, D8_Pearson_r(:,:,3), D8_Pearson_p(:,:,3), ...
    0.05, -1, 1, "Pearson Correlation: NPGO and Max Depth of 8 C Temp")
subplot(2,2,3)
Corr_Map(lat_grid, lon_grid, T814_Spearman_r(:,:,3), T814_Spearman_p(:,:,3), ...
    0.05, -1, 1, "Spearman Correlation: NPGO and Thickness of 8-14 C Temps")
subplot(2,2,4)
Corr_Map(lat_grid, lon_grid, T814_Pearson_r(:,:,3), T814_Pearson_p(:,:,3), ...
    0.05, -1, 1, "Pearson Correlation: NPGO and Thickness of 8-14 C Temps")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Functions
% Function to create maps of mean phase
function Composite_Map(lat_grid_used, lon_grid_used, phase_to_map, min_val, max_val, map_title)
axesm('mercator','MapLatLimit',[10 40],'MapLonLimit',[180 230],'MeridianLabel','on','ParallelLabel','on','Grid','on');
Ig = geoshow(lat_grid_used,lon_grid_used,phase_to_map,'DisplayType','texturemap');
set(Ig,'AlphaData',double(~isnan(phase_to_map)),'AlphaDataMapping','none','FaceAlpha','texturemap'); 
plotm([20 20], [180 230], 'k');
plotm([10 40], [210 210], 'k');
plotm([26 26], [180 210], 'k');
geoshow('landareas.shp','FaceColor',[0.5 0.5 0.5]);
clim([min_val max_val])
cb = colorbar();
% ylabel(cb, sprintf('%s'), variable, 'Rotation', 270); % Matlab does not
% like this
title(sprintf('%s'), map_title)
set(gcf,'renderer','Painters')
tightmap
end

% Function to create map of  max depth of 8 deg C water and thickness of 8-14 deg C
% layer
function Thickness_Map(lat_grid_used, lon_grid_used, depth_to_map, thickness_to_contour, min_val, max_val, map_title)
axesm('mercator','MapLatLimit',[10 40],'MapLonLimit',[180 230],'MeridianLabel','on','ParallelLabel','on','Grid','on');
Ig = geoshow(lat_grid_used,lon_grid_used,depth_to_map,'DisplayType','texturemap');
set(Ig,'AlphaData',double(~isnan(depth_to_map)),'AlphaDataMapping','none','FaceAlpha','texturemap'); 
contourm(lat_grid_used, lon_grid_used, thickness_to_contour, [100 100], 'k', 'LineWidth', 2);
contourm(lat_grid_used, lon_grid_used, thickness_to_contour, 200:100:1000, 'k');
plotm([20 20], [180 230], 'k');
plotm([10 40], [210 210], 'k');
plotm([26 26], [180 210], 'k');
geoshow('landareas.shp','FaceColor',[0.5 0.5 0.5]);
clim([min_val max_val])
cb = colorbar();
% ylabel(cb, sprintf('%s'), variable, 'Rotation', 270); % Matlab does not
% like this
title(sprintf('%s'), map_title)
set(gcf,'renderer','Painters')
tightmap
end

% Function to create maps of the difference from average for composite phases
function Difference_Map(lat_grid_used, lon_grid_used, phase_to_map, min_val, max_val, map_title)
axesm('mercator','MapLatLimit',[10 40],'MapLonLimit',[180 230],'MeridianLabel','on','ParallelLabel','on','Grid','on');
Ig = geoshow(lat_grid_used,lon_grid_used,phase_to_map,'DisplayType','texturemap');
set(Ig,'AlphaData',double(~isnan(phase_to_map)),'AlphaDataMapping','none','FaceAlpha','texturemap'); 
plotm([20 20], [180 230], 'k');
plotm([10 40], [210 210], 'k');
plotm([26 26], [180 210], 'k');
geoshow('landareas.shp','FaceColor',[0.5 0.5 0.5]);
colormap(redblueTecplot)
clim([min_val max_val])
cb = colorbar();
ylabel(cb, 'Difference from Climatology', 'Rotation', 270);
title(sprintf('%s'), map_title)
set(gcf,'renderer','Painters')
tightmap
end

% Function to create map of difference from average for composite phases,
% for depth if 8 deg C and thickness
function Thickness_Diff_Map(lat_grid_used, lon_grid_used, depth_to_map, thickness_to_contour, min_val, max_val, map_title)
axesm('mercator','MapLatLimit',[10 40],'MapLonLimit',[180 230],'MeridianLabel','on','ParallelLabel','on','Grid','on');
Ig = geoshow(lat_grid_used,lon_grid_used,depth_to_map,'DisplayType','texturemap');
set(Ig,'AlphaData',double(~isnan(depth_to_map)),'AlphaDataMapping','none','FaceAlpha','texturemap'); 
contourm(lat_grid_used, lon_grid_used, thickness_to_contour, [0 0], 'k', 'LineWidth', 2);
contourm(lat_grid_used, lon_grid_used, thickness_to_contour, 10:10:100, 'k');
contourm(lat_grid_used, lon_grid_used, thickness_to_contour, -10:-10:-100, 'k:');
plotm([20 20], [180 230], 'k');
plotm([10 40], [210 210], 'k');
plotm([26 26], [180 210], 'k');
geoshow('landareas.shp','FaceColor',[0.5 0.5 0.5]);
colormap(redblueTecplot)
clim([min_val max_val])
cb = colorbar();
% ylabel(cb, sprintf('%s'), variable, 'Rotation', 270); % Matlab does not
% like this
ylabel(cb, 'Difference from Climatology', 'Rotation', 270);
title(sprintf('%s'), map_title)
set(gcf,'renderer','Painters')
tightmap
end

% Function to create map of significant correlations
function Corr_Map(lat_grid_used, lon_grid_used, r_vals, p_vals, sig_level, min_val, max_val, map_title)
r_vals(p_vals >= sig_level) = NaN;
axesm('mercator','MapLatLimit',[10 40],'MapLonLimit',[180 230],'MeridianLabel','on','ParallelLabel','on','Grid','on');
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






























