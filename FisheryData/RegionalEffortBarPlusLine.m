% The purpose of this script is to create a plot that communicates total
% effort over time, plus how that effort is distributed across the fishing
% grounds.  It should be non-confidential

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Load data
Lat = ncread('TotalEffort.nc', 'Latitude');
Lon = ncread('TotalEffort.nc', 'Longitude');
Month = ncread('TotalEffort.nc', 'Month');
Effort = ncread('TotalEffort.nc', 'Total Effort');
Vessels = ncread('TotalVessels.nc', 'Total Vessels');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Wrangle data
% Sum Effort and Vessels over each of the given regions
% Permute so arrays are Lat x Lon x time
Effort = permute(Effort, [2 1 3]);
Vessels = permute(Vessels, [2 1 3]);

% Separate Effort into regions
NW_lat = find(Lat >= 26);
CW_lat = find(Lat >= 20 & Lat < 26);
S_lat = find(Lat >= 10 & Lat < 20);
NE_lat = find(Lat >= 20);
W_lon = find(Lon >= 180 & Lon < 210);
E_lon = find (Lon >= 210 & Lon < 230);

Regional_Effort(1:length(Month), 1:5) = NaN; %Cols are NW, CW, SW, SE, NE
Regional_Vessels(1:length(Month), 1:5) = NaN; %Cols are NW, CW, SW, SE, NE
for m = 1:1:length(Month)
    Regional_Effort(m,1) = sum(Effort(NW_lat, W_lon, m), "all", "omitnan");
    Regional_Effort(m,2) = sum(Effort(CW_lat, W_lon, m), "all", "omitnan");
    Regional_Effort(m,3) = sum(Effort(S_lat, W_lon, m), "all", "omitnan");
    Regional_Effort(m,4) = sum(Effort(S_lat, E_lon, m), "all", "omitnan");
    Regional_Effort(m,5) = sum(Effort(NE_lat, E_lon, m), "all", "omitnan");

    Regional_Vessels(m,1) = sum(Vessels(NW_lat, W_lon, m), "all", "omitnan");
    Regional_Vessels(m,2) = sum(Vessels(CW_lat, W_lon, m), "all", "omitnan");
    Regional_Vessels(m,3) = sum(Vessels(S_lat, W_lon, m), "all", "omitnan");
    Regional_Vessels(m,4) = sum(Vessels(S_lat, E_lon, m), "all", "omitnan");
    Regional_Vessels(m,5) = sum(Vessels(NE_lat, E_lon, m), "all", "omitnan");
end
clear m NW_lat CW_lat S_lat NE_lat W_lon E_lon

% Sum effort and vessels across all regions, too
% Remove singleton dimensions
Total_Effort = squeeze(sum(Effort, [1 2], "omitnan"));
Total_Vessels = squeeze(sum(Vessels, [1 2], "omitnan"));

% Set to zero regional effort when there are less than 3 vessels there that
% month
Regional_Effort(Regional_Vessels < 3) = 0;
Total_Effort(Total_Vessels < 3) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot (these are combined in Illustrator)
figure
bar(Regional_Effort(:,[1 2 3 5])./Total_Effort, 1, 'stacked');
legend('NW', 'CW', 'SW', 'NE');
set(gca,'XTick',1:60:360);
set(gca,'XTickLabel',1995:5:2020);
axis([1 360 0 1]);
pbaspect([2 1 1]);
title('Regional Effort Distribution');

figure
plot(Month, Total_Effort, 'k');
ylabel('Hooks');
set(gca,'XTick',1:60:360);
set(gca,'XTickLabel',1995:5:2020);
axis([1 360 0 6e6]);
pbaspect([2 1 1]);
title('Total Effort');