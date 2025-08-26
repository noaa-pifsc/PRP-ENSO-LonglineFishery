% The purpose of this script is to examine whether there are trends in
% effort over our period of interest, overall and within each region. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Load data
% Matlab needs a "here" package.  Until then, use the path that fits your
% OS.
% Mac:
Lat = ncread('../FisheryData/TotalEffort.nc', 'Latitude');
Lon = ncread('../FisheryData/TotalEffort.nc', 'Longitude');
Month = ncread('../FisheryData/TotalEffort.nc', 'Month');
Effort = ncread('../FisheryData/TotalEffort.nc', 'Total Effort');
% PC:
% Lat = ncread('..\FisheryData\TotalEffort.nc', 'Latitude');
% Lon = ncread('..\FisheryData\TotalEffort.nc', 'Longitude');
% Month = ncread('..\FisheryData\TotalEffort.nc', 'Month');
% Effort = ncread('..\FisheryData\TotalEffort.nc', 'Total Effort');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Wrangle data
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
%%% Linear trends
% Predictor: sequential months, 1–360
x = 1:1:360;

Totlm = fitlm(x, Total_Effort);
NWlm = fitlm(x, Regional_Effort(:,1));
CWlm = fitlm(x, Regional_Effort(:,2));
SWlm = fitlm(x, Regional_Effort(:,3));
SElm = fitlm(x, Regional_Effort(:,4));
NElm = fitlm(x, Regional_Effort(:,5));

% Plot total effort and trend
TotDelta = round((Totlm.Fitted(end) - Totlm.Fitted(1)) / (x(end) - x(1)));
figure
plot(x, Total_Effort, 'k');
hold on
plot(x, Totlm.Fitted,'b');
ylabel('Hooks');
set(gca,'XTick',1:60:360);
set(gca,'XTickLabel',1995:5:2020);
text(5, 5e6, sprintf('+%d Hooks per Month', TotDelta), 'Color', 'b', 'FontSize', 12);
axis([1 360 0 6e6]);
pbaspect([2 1 1]);
title('Total Effort');

% Plot regional effort and trend using fuction at end of script
figure
EffortTrendPlot("Northwest", 1, NWlm, Regional_Effort, x);

figure
EffortTrendPlot("Centralwest", 2, CWlm, Regional_Effort, x);

figure
EffortTrendPlot("Southwest", 3, SWlm, Regional_Effort, x);

figure
EffortTrendPlot("Southeast", 4, SElm, Regional_Effort, x);

figure
EffortTrendPlot("Northeast", 5, NElm, Regional_Effort, x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Anomalies
% There appears to be a seasonal cycle to the data.  We could try removing
% that and plotting the anomaly.
Total_Effort_Climo(1:360,1) = NaN;
Regional_Effort_Climo(1:360,1:5) = NaN;

for m = 1:1:12
    Total_Effort_Climo(m:m:360,1) = mean(Total_Effort(m:m:360,1), "omitnan");

    for r = 1:1:5
        Regional_Effort_Climo(m:m:360,r) = mean(Regional_Effort(m:m:360,r), "omitnan");
    end
end
clear r m

Total_Effort_Anomaly = Total_Effort - Total_Effort_Climo;
Regional_Effort_Anomaly = Regional_Effort - Regional_Effort_Climo;
% This proved not so helpful.  More than anything, it just change the
% values of the min and max.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Functions
% Function to plot monthly effort along with trend and some handy info
function EffortTrendPlot(region_name, region_number, linear_model, fishery, ...
    predictor)
RegDelta = round((linear_model.Fitted(end) - linear_model.Fitted(1)) / ...
    (predictor(end) - predictor(1)));
plot(predictor, fishery(:,region_number), 'k');
hold on
plot(predictor, linear_model.Fitted,'b');
ylabel('Hooks');
set(gca,'XTick',1:60:360);
set(gca,'XTickLabel',1995:5:2020);
text(5, 5e6, sprintf('+%d Hooks per Month', RegDelta), 'Color', 'b', 'FontSize', 12);
axis([1 360 0 6e6]);
pbaspect([2 1 1]);
title(sprintf('%s Region', region_name));
end