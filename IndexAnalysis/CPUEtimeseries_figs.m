% The purpose of this script is to create time series plots of CPUE of our
% four species of interest.  We may also include versions with the climate
% indices included, so I'm doing this from IndexAnalysis.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Load data
% Matlab needs a "here" package.  Until then, use the path that fits your
% OS.
% Mac:
Effort = ncread('../FisheryData/TotalEffort.nc', 'Total Effort');
Bigeye = ncread('../FisheryData/TotalBigeyeCaught.nc', 'Total Bigeye Caught');
Yellowfin = ncread('../FisheryData/TotalYellowfinCaught.nc', 'Total Yellowfin Caught');
Mahi = ncread('../FisheryData/TotalMahiCaught.nc', 'Total Mahi Caught');
Pomfret = ncread('../FisheryData/TotalPomfretCaught.nc', 'Total Pomfret Caught');
ONI = readtable('../ClimateIndices/ONI_withPhases.csv');
PDO = readtable('../ClimateIndices/PDO.csv');
NPGO = readtable('../ClimateIndices/NPGO.csv');
% PC:
% Effort = ncread('..\FisheryData\TotalEffort.nc', 'Total Effort');
% Bigeye = ncread('..\FisheryData\TotalBigeyeCaught.nc', 'Total Bigeye Caught');
% Yellowfin = ncread('..\FisheryData\TotalYellowfinCaught.nc', 'Total Yellowfin Caught');
% Mahi = ncread('..\FisheryData\TotalMahiCaught.nc', 'Total Mahi Caught');
% Pomfret = ncread('..\FisheryData\TotalPomfretCaught.nc', 'Total Pomfret Caught');
% ONI = readtable('..\ClimateIndices\ONI_withPhases.csv');
% PDO = readtable('..\ClimateIndices\PDO.csv');
% NPGO = readtable('..\ClimateIndices\NPGO.csv');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Wrangle data
% Permute so arrays are Lat x Lon x time
Effort = permute(Effort, [2 1 3]);
Bigeye = permute(Bigeye, [2 1 3]);
Yellowfin = permute(Yellowfin, [2 1 3]);
Mahi = permute(Mahi, [2 1 3]);
Pomfret = permute(Pomfret, [2 1 3]);

% Limit climate indices to our period of interest: 1995–2024
ONI_9524 = ONI(ONI.YR >= 1995 & ONI.YR <= 2024,:);
PDO_9524 = PDO(PDO.Year >= 1995 & PDO.Year <= 2024,:);
NPGO_9524 = NPGO(NPGO.YEAR >= 1995 & NPGO.YEAR <= 2024,:);
% Clean up
clear ONI PDO NPGO

% Sum effort and catch across the domain each month
Effort_ts = squeeze(sum(Effort, [1 2], "omitnan"));
Bigeye_ts = squeeze(sum(Bigeye, [1 2], "omitnan"));
Yellowfin_ts = squeeze(sum(Yellowfin, [1 2], "omitnan"));
Pomfret_ts = squeeze(sum(Pomfret, [1 2], "omitnan"));
Mahi_ts = squeeze(sum(Mahi, [1 2], "omitnan"));

% Calculate CPUE
BigeyeCPUE = Bigeye_ts ./ Effort_ts * 1000;
YellowfinCPUE = Yellowfin_ts ./ Effort_ts * 1000;
PomfretCPUE = Pomfret_ts ./ Effort_ts * 1000;
MahiCPUE = Mahi_ts ./ Effort_ts * 1000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CPUE figures
% Plotting all four species together was messy
% So we will make a function to create the plots
figure
CPUE_plot(BigeyeCPUE, 'Bigeye');
% saveas(gcf, 'CPUE_Bigeye.pdf')

figure
CPUE_plot(YellowfinCPUE, 'Yellowfin');
% saveas(gcf, 'CPUE_Yellowfin.pdf')

figure
CPUE_plot(PomfretCPUE, 'Pomfret');
% saveas(gcf, 'CPUE_Pomfret.pdf')

figure
CPUE_plot(MahiCPUE, 'Mahi');
% saveas(gcf, 'CPUE_Mahi.pdf')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Functions
function CPUE_plot(species, fig_title)

plot(species);
set(gca, 'XTick', 1:60:301);
set(gca, 'XTickLabel', 1995:5:2020);
ylabel('CPUE (Fish per 1,000 hooks)');
title(sprintf('%s'), fig_title)
axis([1 360 0 12]);
pbaspect([3 1 1]);

end




