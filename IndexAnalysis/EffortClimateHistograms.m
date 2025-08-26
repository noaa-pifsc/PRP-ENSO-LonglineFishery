% The purpose of this script is to create histograms of climate index value
% with effort.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Load data
% Matlab needs a "here" package.  Until then, use the path that fits your
% OS.
% Mac:
Lat = ncread('../FisheryData/TotalEffort.nc', 'Latitude');
Lon = ncread('../FisheryData/TotalEffort.nc', 'Longitude');
Month = ncread('../FisheryData/TotalEffort.nc', 'Month');
Effort = ncread('../FisheryData/TotalEffort.nc', 'Total Effort');
ONI = readtable('../ClimateIndices/ONI_withPhases.csv');
PDO = readtable('../ClimateIndices/PDO.csv');
NPGO = readtable('../ClimateIndices/NPGO.csv');
% PC:
% Lat = ncread('..\FisheryData\TotalEffort.nc', 'Latitude');
% Lon = ncread('..\FisheryData\TotalEffort.nc', 'Longitude');
% Month = ncread('..\FisheryData\TotalEffort.nc', 'Month');
% Effort = ncread('..\FisheryData\TotalEffort.nc', 'Total Effort');
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

% Sum Effort each of the given regions
% Permute so arrays are Lat x Lon x time
Effort = permute(Effort, [2 1 3]);

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

% Sum effort across all regions, too
% Remove singleton dimensions
Total_Effort = squeeze(sum(Effort, [1 2], "omitnan"));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Histograms
edges = -3.5:0.5:3.5;

% Distribution of climate modes
[N_ONI, edges] = histcounts(ONI_9524.ONI, edges);
[N_PDO, edges] = histcounts(PDO_9524.PDO, edges);
[N_NPGO, edges] = histcounts(NPGO_9524.NPGO, edges);

% Distribution of effort across climate modes, overall and regionally
% Empty array to fill
% region (+total) x bin x mode (ONI, PDO, NPGO)
N_effort(1:6,1:size(edges,2)-1,1:3) = NaN; % region (+total) x bin x mode (ONI, PDO, NPGO)
for g = 1:1:(size(edges,2)-1)
    
    % All but last bin
    if g < (size(edges,2)-1)
        ONI_idx = find(ONI_9524.ONI >= edges(g) & ONI_9524.ONI < edges(g+1));
        PDO_idx = find(PDO_9524.PDO >= edges(g) & PDO_9524.PDO < edges(g+1));
        NPGO_idx = find(NPGO_9524.NPGO >= edges(g) & NPGO_9524.NPGO < edges(g+1));

        % Regional effort
        for r = 1:1:5
            N_effort(r, g, 1) = sum(Regional_Effort(ONI_idx,r), "omitnan");
            N_effort(r, g, 2) = sum(Regional_Effort(PDO_idx,r), "omitnan");
            N_effort(r, g, 3) = sum(Regional_Effort(NPGO_idx,r), "omitnan");
        end

        % Total effort
        N_effort(6, g, 1) = sum(Total_Effort(ONI_idx,1), "omitnan");
        N_effort(6, g, 2) = sum(Total_Effort(PDO_idx,1), "omitnan");
        N_effort(6, g, 3) = sum(Total_Effort(NPGO_idx,1), "omitnan");

        % Last bin
    else
        ONI_idx = find(ONI_9524.ONI >= edges(g) & ONI_9524.ONI <= edges(g+1));
        PDO_idx = find(PDO_9524.PDO >= edges(g) & PDO_9524.PDO <= edges(g+1));
        NPGO_idx = find(NPGO_9524.NPGO >= edges(g) & NPGO_9524.NPGO <= edges(g+1));

        % Regional effort
        for r = 1:1:5
            N_effort(r, g, 1) = sum(Regional_Effort(ONI_idx,r), "omitnan");
            N_effort(r, g, 2) = sum(Regional_Effort(PDO_idx,r), "omitnan");
            N_effort(r, g, 3) = sum(Regional_Effort(NPGO_idx,r), "omitnan");
        end

        % Total effort
        N_effort(6, g, 1) = sum(Total_Effort(ONI_idx,1), "omitnan");
        N_effort(6, g, 2) = sum(Total_Effort(PDO_idx,1), "omitnan");
        N_effort(6, g, 3) = sum(Total_Effort(NPGO_idx,1), "omitnan");
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plots
% Function for this at the end

% ONI
figure
CombinedFigs(N_ONI, 'ONI', 1, N_effort, edges);

% PDO
figure
CombinedFigs(N_PDO, 'PDO', 2,  N_effort, edges);

% NPGO
figure
CombinedFigs(N_NPGO, 'NPGO', 3, N_effort, edges);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Functions
function CombinedFigs(climate_hist, climate_mode, climate_mode_number, effort_hist, bins)
RGB = orderedcolors("gem"); % colors for region (Matlab defaults)
region = ['NW '; 'CW '; 'SW '; 'SE '; 'NE '; 'All']; % for legend
pos_y = 0.325:-0.025:0.2; % for legend
bar(climate_hist./sum(climate_hist), 1); % The 1 is so the bars touch
hold on
for r = 1:1:5
    plot(effort_hist(r,:,climate_mode_number)./sum(effort_hist(r,:,climate_mode_number)), 'Color', RGB(r+1,:), 'LineWidth', 2);
    text(0.5, pos_y(r), region(r,:), 'Color', RGB(r+1,:))
end
plot(effort_hist(6,:,climate_mode_number)./sum(effort_hist(6,:,climate_mode_number)), 'k', 'LineWidth', 2);
text(0.5, pos_y(6), region(6,:), 'Color', 'k')
set(gca, 'XTick', 0.5:1:14.5);
set(gca, 'XTickLabel', bins);
xlabel(sprintf('%s', climate_mode));
ylabel('Proportion of Months (bars) and Effort (lines)');
ylim([0 0.35]);
grid on
pbaspect([2 1 1]);
end



































