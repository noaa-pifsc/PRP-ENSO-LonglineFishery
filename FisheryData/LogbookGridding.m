tic
% The purpose of this script is to grid logbook data.  It reads in
% set-level logbook data (see LogbookAccess.Rmd and LogbookCombine.Rmd) and
% prodces data gridded at 1-deg x 1-deg and monthly resolution.  The grid
% matches that used in OceanDataRegrid.jnl.  This means that all our data
% are on the same grid.
%
% Approach to gridding:
% Grid cells are centered on half-degrees and grid edges are at whole
% degrees.  Grids are anchored in the lower left or southwest corner
% (i.e., the >= and < operators are used).
%
% Data fields:
%   Effort: Total number of hooks set in a month
%   Sets: Total number of sets in a month
%   Bigeye: Total number of bigeye kept + discarded
%   Swordfish: Total number of swordfish kept + discarded
%   Mahi: Total number of mahi kept + discarded
%   Yellowfin: Total number of yellowfin kept + discarded
%   Pomfret: Total number of pomfret kept + discarded
%   Vessels: Number of vessels that set effort in a given cell, aids in 
%            determining whether results can be shared publicly.  Unique
%            permit numbers are assumed to represent unique vessels.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Load data
% Catchability (to ensure grids match)
% Matlab needs a "here" package.  Until then, use the path that fits your
% OS.
% Mac:
lat = ncread('../OceanData/BigeyeCatchability.nc', 'Latitude');
lon = ncread('../OceanData/BigeyeCatchability.nc', 'Longitude');
month = ncread('../OceanData/BigeyeCatchability.nc', 'Month');
% PC:
% lat = ncread('..\OceanData\BigeyeCatchability.nc', 'Latitude');
% lon = ncread('..\OceanData\BigeyeCatchability.nc', 'Longitude');
% month = ncread('..\OceanData\BigeyeCatchability.nc', 'Month');

% Set-level data
DeepSets = readtable('DeepSets_AllSpecies.csv');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Data wrangling
% The logbook longitudes are in positive degrees west (i.e., 150 = 150W)
% There are no longitudes < 0 or >= 180
% Converting to 360-longitude
DeepSets.BH_LON = 360 - DeepSets.BH_LON;

% Empty arrays to fill
Effort = NaN(length(lat), length(lon), length(month));
Sets = NaN(length(lat), length(lon), length(month));
Bigeye = NaN(length(lat), length(lon), length(month));
Swordfish = NaN(length(lat), length(lon), length(month));
Mahi = NaN(length(lat), length(lon), length(month));
Yellowfin = NaN(length(lat), length(lon), length(month));
Pomfret = NaN(length(lat), length(lon), length(month));
Vessels = NaN(length(lat), length(lon), length(month));
PMUS = NaN(length(lat), length(lon), length(month));
Richness = NaN(length(lat), length(lon), length(month));
Shannon = NaN(length(lat), length(lon), length(month));
Simpson = NaN(length(lat), length(lon), length(month));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Gridding
% Loop through each year, month, lat, and lon
m = 1; % counter
for yr = 1995:1:2024
    for mo = 1:1:12
        for r = 1:1:length(lat)
            for c = 1:1:length(lon)
                idx = find(DeepSets.BH_YR == yr & ...
                           DeepSets.BH_MON == mo & ...
                           DeepSets.BH_LAT >= (lat(r) - 0.5) & ...
                           DeepSets.BH_LAT < (lat(r) + 0.5) & ...
                           DeepSets.BH_LON >= (lon(c) - 0.5) & ...
                           DeepSets.BH_LON < (lon(c) + 0.5));

                Effort(r,c,m) = sum(DeepSets.HOOKSSET(idx), "omitnan");
                Sets(r,c,m) = length(idx); % each row is a set
                Bigeye(r,c,m) = sum(DeepSets.NUMKEPT_16(idx), "omitnan") + ...
                                sum(DeepSets.NUMRELEASED_16(idx), "omitnan");
                Swordfish(r,c,m) = sum(DeepSets.NUMKEPT_6(idx), "omitnan") + ...
                                   sum(DeepSets.NUMRELEASED_6(idx), "omitnan");
                Mahi(r,c,m) = sum(DeepSets.NUMKEPT_11(idx), "omitnan") + ...
                              sum(DeepSets.NUMRELEASED_11(idx), "omitnan");
                Yellowfin(r,c,m) = sum(DeepSets.NUMKEPT_17(idx), "omitnan") + ...
                                   sum(DeepSets.NUMRELEASED_17(idx), "omitnan");
                Pomfret(r,c,m) = sum(DeepSets.NUMKEPT_21(idx), "omitnan") + ...
                                 sum(DeepSets.NUMRELEASED_21(idx), "omitnan");
                Vessels(r,c,m) = length(unique(DeepSets.PERMITNUM(idx)));

                % For the PMUS species and the diversity indices, we're
                % going to make a smaller table that we can operate over
                SmallerTable = DeepSets(idx,{'NUMKEPT_1','NUMRELEASED_1', ...
                                             'NUMKEPT_2','NUMRELEASED_2', ...
                                             'NUMKEPT_3','NUMRELEASED_3', ...
                                             'NUMKEPT_4','NUMRELEASED_4', ...
                                             'NUMKEPT_5','NUMRELEASED_5', ...
                                             'NUMKEPT_6','NUMRELEASED_6', ...
                                             'NUMFINNED_7','NUMKEPT_7','NUMRELEASED_7', ...
                                             'NUMFINNED_8','NUMKEPT_8','NUMRELEASED_8', ...
                                             'NUMFINNED_9','NUMKEPT_9','NUMRELEASED_9', ...
                                             'NUMKEPT_11','NUMRELEASED_11', ...
                                             'NUMKEPT_12','NUMRELEASED_12', ...
                                             'NUMKEPT_13','NUMRELEASED_13', ...
                                             'NUMKEPT_15','NUMRELEASED_15', ...
                                             'NUMKEPT_16','NUMRELEASED_16', ...
                                             'NUMKEPT_17','NUMRELEASED_17', ...
                                             'NUMKEPT_18','NUMRELEASED_18', ...
                                             'NUMKEPT_19','NUMRELEASED_19', ...
                                             'NUMKEPT_20','NUMRELEASED_20', ...
                                             'NUMKEPT_21','NUMRELEASED_21', ...
                                             'NUMKEPT_22','NUMRELEASED_22', ...
                                             'NUMFINNED_24','NUMKEPT_24','NUMRELEASED_24', ...
                                             'NUMFINNED_25','NUMKEPT_25','NUMRELEASED_25'});
                
                % Sum SmallerTable columns
                ST_colsums = sum(SmallerTable, 1, "omitnan");

                % Sum across all columns
                ST_sum = sum(SmallerTable, 'all');
                
                % PMUS total catch
                PMUS(r,c,m) = ST_sum{1,1};

                % Richness
                % This is tedious, but I can't think of a better way to
                % ensure we don't double- or triple-count species
                % Suggestions welcome!
                SP_1 = 0; SP_2 = 0; SP_3 = 0; SP_4 = 0; SP_5 = 0; SP_6 = 0;
                SP_7 = 0; SP_8 = 0; SP_9 = 0; SP_11 = 0; SP_12 = 0; SP_13 = 0;
                SP_15 = 0; SP_16 = 0; SP_17 = 0; SP_18 = 0; SP_19 = 0;
                SP_20 = 0; SP_21 = 0; SP_22 = 0; SP_24 = 0; SP_25 = 0;

                if ST_colsums.NUMKEPT_1 > 0 || ST_colsums.NUMRELEASED_1 > 0
                    SP_1 = 1;
                end
                if ST_colsums.NUMKEPT_2 > 0 || ST_colsums.NUMRELEASED_2 > 0
                    SP_2 = 1;
                end
                if ST_colsums.NUMKEPT_3 > 0 || ST_colsums.NUMRELEASED_3 > 0
                    SP_3 = 1;
                end
                if ST_colsums.NUMKEPT_4 > 0 || ST_colsums.NUMRELEASED_4 > 0
                    SP_4 = 1;
                end
                if ST_colsums.NUMKEPT_5 > 0 || ST_colsums.NUMRELEASED_5 > 0
                    SP_5 = 1;
                end
                if ST_colsums.NUMKEPT_6 > 0 || ST_colsums.NUMRELEASED_6 > 0
                    SP_6 = 1;
                end
                if ST_colsums.NUMFINNED_7 > 0 || ST_colsums.NUMKEPT_7 > 0 || ST_colsums.NUMRELEASED_7 > 0
                    SP_7 = 1;
                end
                if ST_colsums.NUMFINNED_8 > 0 || ST_colsums.NUMKEPT_8 > 0 || ST_colsums.NUMRELEASED_8 > 0
                    SP_8 = 1;
                end
                if ST_colsums.NUMFINNED_9 > 0 || ST_colsums.NUMKEPT_9 > 0 || ST_colsums.NUMRELEASED_9 > 0
                    SP_9 = 1;
                end
                if ST_colsums.NUMKEPT_11 > 0 || ST_colsums.NUMRELEASED_11 > 0
                    SP_11 = 1;
                end
                if ST_colsums.NUMKEPT_12 > 0 || ST_colsums.NUMRELEASED_12 > 0
                    SP_12 = 1;
                end
                if ST_colsums.NUMKEPT_13 > 0 || ST_colsums.NUMRELEASED_13 > 0
                    SP_13 = 1;
                end
                if ST_colsums.NUMKEPT_15 > 0 || ST_colsums.NUMRELEASED_15 > 0
                    SP_15 = 1;
                end
                if ST_colsums.NUMKEPT_16 > 0 || ST_colsums.NUMRELEASED_16 > 0
                    SP_16 = 1;
                end
                if ST_colsums.NUMKEPT_17 > 0 || ST_colsums.NUMRELEASED_17 > 0
                    SP_17 = 1;
                end
                if ST_colsums.NUMKEPT_18 > 0 || ST_colsums.NUMRELEASED_18 > 0
                    SP_18 = 1;
                end
                if ST_colsums.NUMKEPT_19 > 0 || ST_colsums.NUMRELEASED_19 > 0
                    SP_19 = 1;
                end
                if ST_colsums.NUMKEPT_20 > 0 || ST_colsums.NUMRELEASED_20 > 0
                    SP_20 = 1;
                end
                if ST_colsums.NUMKEPT_21 > 0 || ST_colsums.NUMRELEASED_21 > 0
                    SP_21 = 1;
                end
                if ST_colsums.NUMKEPT_22 > 0 || ST_colsums.NUMRELEASED_22 > 0
                    SP_22 = 1;
                end
                if ST_colsums.NUMFINNED_24 > 0 || ST_colsums.NUMKEPT_24 > 0 || ST_colsums.NUMRELEASED_24 > 0
                    SP_24 = 1;
                end
                if ST_colsums.NUMFINNED_25 > 0 || ST_colsums.NUMKEPT_25 > 0 || ST_colsums.NUMRELEASED_25 > 0
                    SP_25 = 1;
                end

                Richness(r,c,m) = SP_1 + SP_2 + SP_3 + SP_4 + SP_5 + SP_6 + ...
                                  SP_7 + SP_8 + SP_9 + SP_11 + SP_12 + ...
                                  SP_13 + SP_15 + SP_16 + SP_17 + SP_18 + ...
                                  SP_19 + SP_20 + SP_21 + SP_22 + SP_24 + ...
                                  SP_25;

                % Shannon
                GT0 = ST_colsums > 0; % Removing zeros to keep the math clean
                ST_PropCatch = ST_colsums(1,GT0{1,:})./ST_sum{1,1}; 
                LN_PropCatch = log(ST_PropCatch);
                H = sum(ST_PropCatch.*LN_PropCatch, 'all');
                if sum(size(H)) > 0
                    Shannon(r,c,m) = H{1,1}*-1;
                end

                % Simpson
                D = sum(ST_PropCatch.^2, 'all');
                if sum(size(D)) > 0
                    Simpson(r,c,m) = 1-D{1,1};
                end

                clear idx SmallerTable ST* SP*

            end
        end
        m = m + 1; % update counter

    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% More data wrangling
% Permute, so that data are lon x lat x time
Effort = permute(Effort, [2 1 3]);
Sets = permute(Sets, [2 1 3]);
Bigeye = permute(Bigeye, [2 1 3]);
Swordfish = permute(Swordfish, [2 1 3]);
Mahi = permute(Mahi, [2 1 3]);
Yellowfin = permute(Yellowfin, [2 1 3]);
Pomfret = permute(Pomfret, [2 1 3]);
Vessels = permute(Vessels, [2 1 3]);
PMUS = permute(PMUS, [2 1 3]);
Richness = permute(Richness, [2 1 3]);
Shannon = permute(Shannon, [2 1 3]);
Simpson = permute(Simpson, [2 1 3]);

% Create lat and lon grids
% Horizontally stack column with latitudes
lat_grid = repmat(lat,1,size(lon,1));
% Transpose column with longitudes to a row and vertically stack
lon_grid = repmat(lon',size(lat,1),1);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Save data
% This approach seems overly tedious because common dimensions and 
% variables have to be recreated for each file...
% And we can only do one file at a time.

% Create file, noting NOT to overwrite an exsiting file
ncid1 = netcdf.create('TotalEffort.nc','NOCLOBBER');

% Define dimensions
dimid_lon1 = netcdf.defDim(ncid1,'Longitude',size(lon_grid,2));
dimid_lat1 = netcdf.defDim(ncid1,'Latitude',size(lat_grid,1));
dimid_mon1 = netcdf.defDim(ncid1,'Month',size(Effort,3)); 

% Define variables
varid_lon1 = netcdf.defVar(ncid1,'Longitude','double',dimid_lon1);
varid_lat1 = netcdf.defVar(ncid1,'Latitude','double',dimid_lat1);
varid_mon1 = netcdf.defVar(ncid1,'Month','double',dimid_mon1);

varid_DS = netcdf.defVar(ncid1,'Total Effort','double',[dimid_lon1 dimid_lat1 dimid_mon1]);

% Define attributes
netcdf.putAtt(ncid1,varid_lon1,'standard_name','Longitude');
netcdf.putAtt(ncid1,varid_lon1,'units','Degrees East');
netcdf.putAtt(ncid1,varid_lon1,'reference','Grid Center');

netcdf.putAtt(ncid1,varid_lat1,'standard_name','Latitude');
netcdf.putAtt(ncid1,varid_lat1,'units','Degrees North');
netcdf.putAtt(ncid1,varid_lat1,'reference','Grid Center');

netcdf.putAtt(ncid1,varid_mon1,'standard_name','Month');
netcdf.putAtt(ncid1,varid_mon1,'units','Months Since Dec 1994');

netcdf.putAtt(ncid1,varid_DS,'standard_name','Total Effort');
netcdf.putAtt(ncid1,varid_DS,'units','Hooks Set');

netcdf.endDef(ncid1)
% netcdf.reDef(ncid1) in case it's necessary to reenter define mode.

% Put the data in the file
netcdf.putVar(ncid1,varid_lon1,lon);
netcdf.putVar(ncid1,varid_lat1,lat);
netcdf.putVar(ncid1,varid_mon1,month); % time since Dec 1994 not 1992
netcdf.putVar(ncid1,varid_DS,Effort);

% Close the files so they can be used
netcdf.close(ncid1)

% Clean up, so we can recycle code (sigh...)
clear ncid* dimid* varid*

%%% Create the other files
% Lines that change across files are indented 
    ncid1 = netcdf.create('TotalSets.nc','NOCLOBBER');
dimid_lon1 = netcdf.defDim(ncid1,'Longitude',size(lon_grid,2));
dimid_lat1 = netcdf.defDim(ncid1,'Latitude',size(lat_grid,1));
dimid_mon1 = netcdf.defDim(ncid1,'Month',size(Effort,3));
varid_lon1 = netcdf.defVar(ncid1,'Longitude','double',dimid_lon1);
varid_lat1 = netcdf.defVar(ncid1,'Latitude','double',dimid_lat1);
varid_mon1 = netcdf.defVar(ncid1,'Month','double',dimid_mon1);
    varid_DS = netcdf.defVar(ncid1,'Total Sets','double',[dimid_lon1 dimid_lat1 dimid_mon1]);
netcdf.putAtt(ncid1,varid_lon1,'standard_name','Longitude');
netcdf.putAtt(ncid1,varid_lon1,'units','Degrees East');
netcdf.putAtt(ncid1,varid_lon1,'reference','Grid Center');
netcdf.putAtt(ncid1,varid_lat1,'standard_name','Latitude');
netcdf.putAtt(ncid1,varid_lat1,'units','Degrees North');
netcdf.putAtt(ncid1,varid_lat1,'reference','Grid Center');
netcdf.putAtt(ncid1,varid_mon1,'standard_name','Month');
netcdf.putAtt(ncid1,varid_mon1,'units','Months Since Dec 1994');
    netcdf.putAtt(ncid1,varid_DS,'standard_name','Total Sets');
    netcdf.putAtt(ncid1,varid_DS,'units','Sets');
netcdf.endDef(ncid1)
netcdf.putVar(ncid1,varid_lon1,lon);
netcdf.putVar(ncid1,varid_lat1,lat);
netcdf.putVar(ncid1,varid_mon1,month); 
    netcdf.putVar(ncid1,varid_DS,Sets);
netcdf.close(ncid1)
clear ncid* dimid* varid*

    ncid1 = netcdf.create('TotalBigeyeCaught.nc','NOCLOBBER');
dimid_lon1 = netcdf.defDim(ncid1,'Longitude',size(lon_grid,2));
dimid_lat1 = netcdf.defDim(ncid1,'Latitude',size(lat_grid,1));
dimid_mon1 = netcdf.defDim(ncid1,'Month',size(Effort,3));
varid_lon1 = netcdf.defVar(ncid1,'Longitude','double',dimid_lon1);
varid_lat1 = netcdf.defVar(ncid1,'Latitude','double',dimid_lat1);
varid_mon1 = netcdf.defVar(ncid1,'Month','double',dimid_mon1);
    varid_DS = netcdf.defVar(ncid1,'Total Bigeye Caught','double',[dimid_lon1 dimid_lat1 dimid_mon1]);
netcdf.putAtt(ncid1,varid_lon1,'standard_name','Longitude');
netcdf.putAtt(ncid1,varid_lon1,'units','Degrees East');
netcdf.putAtt(ncid1,varid_lon1,'reference','Grid Center');
netcdf.putAtt(ncid1,varid_lat1,'standard_name','Latitude');
netcdf.putAtt(ncid1,varid_lat1,'units','Degrees North');
netcdf.putAtt(ncid1,varid_lat1,'reference','Grid Center');
netcdf.putAtt(ncid1,varid_mon1,'standard_name','Month');
netcdf.putAtt(ncid1,varid_mon1,'units','Months Since Dec 1994');
    netcdf.putAtt(ncid1,varid_DS,'standard_name','Total Bigeye Caught');
    netcdf.putAtt(ncid1,varid_DS,'units','Fish');
netcdf.endDef(ncid1)
netcdf.putVar(ncid1,varid_lon1,lon);
netcdf.putVar(ncid1,varid_lat1,lat);
netcdf.putVar(ncid1,varid_mon1,month); 
    netcdf.putVar(ncid1,varid_DS,Bigeye);
netcdf.close(ncid1)
clear ncid* dimid* varid*

    ncid1 = netcdf.create('TotalSwordfishCaught.nc','NOCLOBBER');
dimid_lon1 = netcdf.defDim(ncid1,'Longitude',size(lon_grid,2));
dimid_lat1 = netcdf.defDim(ncid1,'Latitude',size(lat_grid,1));
dimid_mon1 = netcdf.defDim(ncid1,'Month',size(Effort,3));
varid_lon1 = netcdf.defVar(ncid1,'Longitude','double',dimid_lon1);
varid_lat1 = netcdf.defVar(ncid1,'Latitude','double',dimid_lat1);
varid_mon1 = netcdf.defVar(ncid1,'Month','double',dimid_mon1);
    varid_DS = netcdf.defVar(ncid1,'Total Swordfish Caught','double',[dimid_lon1 dimid_lat1 dimid_mon1]);
netcdf.putAtt(ncid1,varid_lon1,'standard_name','Longitude');
netcdf.putAtt(ncid1,varid_lon1,'units','Degrees East');
netcdf.putAtt(ncid1,varid_lon1,'reference','Grid Center');
netcdf.putAtt(ncid1,varid_lat1,'standard_name','Latitude');
netcdf.putAtt(ncid1,varid_lat1,'units','Degrees North');
netcdf.putAtt(ncid1,varid_lat1,'reference','Grid Center');
netcdf.putAtt(ncid1,varid_mon1,'standard_name','Month');
netcdf.putAtt(ncid1,varid_mon1,'units','Months Since Dec 1994');
    netcdf.putAtt(ncid1,varid_DS,'standard_name','Total Swordfish Caught');
    netcdf.putAtt(ncid1,varid_DS,'units','Fish');
netcdf.endDef(ncid1)
netcdf.putVar(ncid1,varid_lon1,lon);
netcdf.putVar(ncid1,varid_lat1,lat);
netcdf.putVar(ncid1,varid_mon1,month); 
    netcdf.putVar(ncid1,varid_DS,Swordfish);
netcdf.close(ncid1)
clear ncid* dimid* varid*

    ncid1 = netcdf.create('TotalMahiCaught.nc','NOCLOBBER');
dimid_lon1 = netcdf.defDim(ncid1,'Longitude',size(lon_grid,2));
dimid_lat1 = netcdf.defDim(ncid1,'Latitude',size(lat_grid,1));
dimid_mon1 = netcdf.defDim(ncid1,'Month',size(Effort,3));
varid_lon1 = netcdf.defVar(ncid1,'Longitude','double',dimid_lon1);
varid_lat1 = netcdf.defVar(ncid1,'Latitude','double',dimid_lat1);
varid_mon1 = netcdf.defVar(ncid1,'Month','double',dimid_mon1);
    varid_DS = netcdf.defVar(ncid1,'Total Mahi Caught','double',[dimid_lon1 dimid_lat1 dimid_mon1]);
netcdf.putAtt(ncid1,varid_lon1,'standard_name','Longitude');
netcdf.putAtt(ncid1,varid_lon1,'units','Degrees East');
netcdf.putAtt(ncid1,varid_lon1,'reference','Grid Center');
netcdf.putAtt(ncid1,varid_lat1,'standard_name','Latitude');
netcdf.putAtt(ncid1,varid_lat1,'units','Degrees North');
netcdf.putAtt(ncid1,varid_lat1,'reference','Grid Center');
netcdf.putAtt(ncid1,varid_mon1,'standard_name','Month');
netcdf.putAtt(ncid1,varid_mon1,'units','Months Since Dec 1994');
    netcdf.putAtt(ncid1,varid_DS,'standard_name','Total Mahi Caught');
    netcdf.putAtt(ncid1,varid_DS,'units','Fish');
netcdf.endDef(ncid1)
netcdf.putVar(ncid1,varid_lon1,lon);
netcdf.putVar(ncid1,varid_lat1,lat);
netcdf.putVar(ncid1,varid_mon1,month); 
    netcdf.putVar(ncid1,varid_DS,Mahi);
netcdf.close(ncid1)
clear ncid* dimid* varid*

    ncid1 = netcdf.create('TotalYellowfinCaught.nc','NOCLOBBER');
dimid_lon1 = netcdf.defDim(ncid1,'Longitude',size(lon_grid,2));
dimid_lat1 = netcdf.defDim(ncid1,'Latitude',size(lat_grid,1));
dimid_mon1 = netcdf.defDim(ncid1,'Month',size(Effort,3));
varid_lon1 = netcdf.defVar(ncid1,'Longitude','double',dimid_lon1);
varid_lat1 = netcdf.defVar(ncid1,'Latitude','double',dimid_lat1);
varid_mon1 = netcdf.defVar(ncid1,'Month','double',dimid_mon1);
    varid_DS = netcdf.defVar(ncid1,'Total Yellowfin Caught','double',[dimid_lon1 dimid_lat1 dimid_mon1]);
netcdf.putAtt(ncid1,varid_lon1,'standard_name','Longitude');
netcdf.putAtt(ncid1,varid_lon1,'units','Degrees East');
netcdf.putAtt(ncid1,varid_lon1,'reference','Grid Center');
netcdf.putAtt(ncid1,varid_lat1,'standard_name','Latitude');
netcdf.putAtt(ncid1,varid_lat1,'units','Degrees North');
netcdf.putAtt(ncid1,varid_lat1,'reference','Grid Center');
netcdf.putAtt(ncid1,varid_mon1,'standard_name','Month');
netcdf.putAtt(ncid1,varid_mon1,'units','Months Since Dec 1994');
    netcdf.putAtt(ncid1,varid_DS,'standard_name','Total Yellowfin Caught');
    netcdf.putAtt(ncid1,varid_DS,'units','Fish');
netcdf.endDef(ncid1)
netcdf.putVar(ncid1,varid_lon1,lon);
netcdf.putVar(ncid1,varid_lat1,lat);
netcdf.putVar(ncid1,varid_mon1,month); 
    netcdf.putVar(ncid1,varid_DS,Yellowfin);
netcdf.close(ncid1)
clear ncid* dimid* varid*

    ncid1 = netcdf.create('TotalPomfretCaught.nc','NOCLOBBER');
dimid_lon1 = netcdf.defDim(ncid1,'Longitude',size(lon_grid,2));
dimid_lat1 = netcdf.defDim(ncid1,'Latitude',size(lat_grid,1));
dimid_mon1 = netcdf.defDim(ncid1,'Month',size(Effort,3));
varid_lon1 = netcdf.defVar(ncid1,'Longitude','double',dimid_lon1);
varid_lat1 = netcdf.defVar(ncid1,'Latitude','double',dimid_lat1);
varid_mon1 = netcdf.defVar(ncid1,'Month','double',dimid_mon1);
    varid_DS = netcdf.defVar(ncid1,'Total Pomfret Caught','double',[dimid_lon1 dimid_lat1 dimid_mon1]);
netcdf.putAtt(ncid1,varid_lon1,'standard_name','Longitude');
netcdf.putAtt(ncid1,varid_lon1,'units','Degrees East');
netcdf.putAtt(ncid1,varid_lon1,'reference','Grid Center');
netcdf.putAtt(ncid1,varid_lat1,'standard_name','Latitude');
netcdf.putAtt(ncid1,varid_lat1,'units','Degrees North');
netcdf.putAtt(ncid1,varid_lat1,'reference','Grid Center');
netcdf.putAtt(ncid1,varid_mon1,'standard_name','Month');
netcdf.putAtt(ncid1,varid_mon1,'units','Months Since Dec 1994');
    netcdf.putAtt(ncid1,varid_DS,'standard_name','Total Pomfret Caught');
    netcdf.putAtt(ncid1,varid_DS,'units','Fish');
netcdf.endDef(ncid1)
netcdf.putVar(ncid1,varid_lon1,lon);
netcdf.putVar(ncid1,varid_lat1,lat);
netcdf.putVar(ncid1,varid_mon1,month); 
    netcdf.putVar(ncid1,varid_DS,Pomfret);
netcdf.close(ncid1)
clear ncid* dimid* varid*

    ncid1 = netcdf.create('TotalVessels.nc','NOCLOBBER');
dimid_lon1 = netcdf.defDim(ncid1,'Longitude',size(lon_grid,2));
dimid_lat1 = netcdf.defDim(ncid1,'Latitude',size(lat_grid,1));
dimid_mon1 = netcdf.defDim(ncid1,'Month',size(Effort,3));
varid_lon1 = netcdf.defVar(ncid1,'Longitude','double',dimid_lon1);
varid_lat1 = netcdf.defVar(ncid1,'Latitude','double',dimid_lat1);
varid_mon1 = netcdf.defVar(ncid1,'Month','double',dimid_mon1);
    varid_DS = netcdf.defVar(ncid1,'Total Vessels','double',[dimid_lon1 dimid_lat1 dimid_mon1]);
netcdf.putAtt(ncid1,varid_lon1,'standard_name','Longitude');
netcdf.putAtt(ncid1,varid_lon1,'units','Degrees East');
netcdf.putAtt(ncid1,varid_lon1,'reference','Grid Center');
netcdf.putAtt(ncid1,varid_lat1,'standard_name','Latitude');
netcdf.putAtt(ncid1,varid_lat1,'units','Degrees North');
netcdf.putAtt(ncid1,varid_lat1,'reference','Grid Center');
netcdf.putAtt(ncid1,varid_mon1,'standard_name','Month');
netcdf.putAtt(ncid1,varid_mon1,'units','Months Since Dec 1994');
    netcdf.putAtt(ncid1,varid_DS,'standard_name','Total Vessels');
    netcdf.putAtt(ncid1,varid_DS,'units','Vessels');
netcdf.endDef(ncid1)
netcdf.putVar(ncid1,varid_lon1,lon);
netcdf.putVar(ncid1,varid_lat1,lat);
netcdf.putVar(ncid1,varid_mon1,month); 
    netcdf.putVar(ncid1,varid_DS,Vessels);
netcdf.close(ncid1)
clear ncid* dimid* varid*

    ncid1 = netcdf.create('TotalPMUS.nc','NOCLOBBER');
dimid_lon1 = netcdf.defDim(ncid1,'Longitude',size(lon_grid,2));
dimid_lat1 = netcdf.defDim(ncid1,'Latitude',size(lat_grid,1));
dimid_mon1 = netcdf.defDim(ncid1,'Month',size(Effort,3));
varid_lon1 = netcdf.defVar(ncid1,'Longitude','double',dimid_lon1);
varid_lat1 = netcdf.defVar(ncid1,'Latitude','double',dimid_lat1);
varid_mon1 = netcdf.defVar(ncid1,'Month','double',dimid_mon1);
    varid_DS = netcdf.defVar(ncid1,'Total PMUS Caught','double',[dimid_lon1 dimid_lat1 dimid_mon1]);
netcdf.putAtt(ncid1,varid_lon1,'standard_name','Longitude');
netcdf.putAtt(ncid1,varid_lon1,'units','Degrees East');
netcdf.putAtt(ncid1,varid_lon1,'reference','Grid Center');
netcdf.putAtt(ncid1,varid_lat1,'standard_name','Latitude');
netcdf.putAtt(ncid1,varid_lat1,'units','Degrees North');
netcdf.putAtt(ncid1,varid_lat1,'reference','Grid Center');
netcdf.putAtt(ncid1,varid_mon1,'standard_name','Month');
netcdf.putAtt(ncid1,varid_mon1,'units','Months Since Dec 1994');
    netcdf.putAtt(ncid1,varid_DS,'standard_name','Total PMUS Caught');
    netcdf.putAtt(ncid1,varid_DS,'units','Fish');
netcdf.endDef(ncid1)
netcdf.putVar(ncid1,varid_lon1,lon);
netcdf.putVar(ncid1,varid_lat1,lat);
netcdf.putVar(ncid1,varid_mon1,month); 
    netcdf.putVar(ncid1,varid_DS,PMUS);
netcdf.close(ncid1)
clear ncid* dimid* varid*

    ncid1 = netcdf.create('SpeciesRichness.nc','NOCLOBBER');
dimid_lon1 = netcdf.defDim(ncid1,'Longitude',size(lon_grid,2));
dimid_lat1 = netcdf.defDim(ncid1,'Latitude',size(lat_grid,1));
dimid_mon1 = netcdf.defDim(ncid1,'Month',size(Effort,3));
varid_lon1 = netcdf.defVar(ncid1,'Longitude','double',dimid_lon1);
varid_lat1 = netcdf.defVar(ncid1,'Latitude','double',dimid_lat1);
varid_mon1 = netcdf.defVar(ncid1,'Month','double',dimid_mon1);
    varid_DS = netcdf.defVar(ncid1,'Species Richness','double',[dimid_lon1 dimid_lat1 dimid_mon1]);
netcdf.putAtt(ncid1,varid_lon1,'standard_name','Longitude');
netcdf.putAtt(ncid1,varid_lon1,'units','Degrees East');
netcdf.putAtt(ncid1,varid_lon1,'reference','Grid Center');
netcdf.putAtt(ncid1,varid_lat1,'standard_name','Latitude');
netcdf.putAtt(ncid1,varid_lat1,'units','Degrees North');
netcdf.putAtt(ncid1,varid_lat1,'reference','Grid Center');
netcdf.putAtt(ncid1,varid_mon1,'standard_name','Month');
netcdf.putAtt(ncid1,varid_mon1,'units','Months Since Dec 1994');
    netcdf.putAtt(ncid1,varid_DS,'standard_name','Species Richness');
    netcdf.putAtt(ncid1,varid_DS,'units','Species');
netcdf.endDef(ncid1)
netcdf.putVar(ncid1,varid_lon1,lon);
netcdf.putVar(ncid1,varid_lat1,lat);
netcdf.putVar(ncid1,varid_mon1,month); 
    netcdf.putVar(ncid1,varid_DS,Richness);
netcdf.close(ncid1)
clear ncid* dimid* varid*

    ncid1 = netcdf.create('ShannonIndex.nc','NOCLOBBER');
dimid_lon1 = netcdf.defDim(ncid1,'Longitude',size(lon_grid,2));
dimid_lat1 = netcdf.defDim(ncid1,'Latitude',size(lat_grid,1));
dimid_mon1 = netcdf.defDim(ncid1,'Month',size(Effort,3));
varid_lon1 = netcdf.defVar(ncid1,'Longitude','double',dimid_lon1);
varid_lat1 = netcdf.defVar(ncid1,'Latitude','double',dimid_lat1);
varid_mon1 = netcdf.defVar(ncid1,'Month','double',dimid_mon1);
    varid_DS = netcdf.defVar(ncid1,'Shannon Index','double',[dimid_lon1 dimid_lat1 dimid_mon1]);
netcdf.putAtt(ncid1,varid_lon1,'standard_name','Longitude');
netcdf.putAtt(ncid1,varid_lon1,'units','Degrees East');
netcdf.putAtt(ncid1,varid_lon1,'reference','Grid Center');
netcdf.putAtt(ncid1,varid_lat1,'standard_name','Latitude');
netcdf.putAtt(ncid1,varid_lat1,'units','Degrees North');
netcdf.putAtt(ncid1,varid_lat1,'reference','Grid Center');
netcdf.putAtt(ncid1,varid_mon1,'standard_name','Month');
netcdf.putAtt(ncid1,varid_mon1,'units','Months Since Dec 1994');
    netcdf.putAtt(ncid1,varid_DS,'standard_name','Shannon Index');
    netcdf.putAtt(ncid1,varid_DS,'units','Unitless');
netcdf.endDef(ncid1)
netcdf.putVar(ncid1,varid_lon1,lon);
netcdf.putVar(ncid1,varid_lat1,lat);
netcdf.putVar(ncid1,varid_mon1,month); 
    netcdf.putVar(ncid1,varid_DS,Shannon);
netcdf.close(ncid1)
clear ncid* dimid* varid*

    ncid1 = netcdf.create('SimpsonIndex.nc','NOCLOBBER');
dimid_lon1 = netcdf.defDim(ncid1,'Longitude',size(lon_grid,2));
dimid_lat1 = netcdf.defDim(ncid1,'Latitude',size(lat_grid,1));
dimid_mon1 = netcdf.defDim(ncid1,'Month',size(Effort,3));
varid_lon1 = netcdf.defVar(ncid1,'Longitude','double',dimid_lon1);
varid_lat1 = netcdf.defVar(ncid1,'Latitude','double',dimid_lat1);
varid_mon1 = netcdf.defVar(ncid1,'Month','double',dimid_mon1);
    varid_DS = netcdf.defVar(ncid1,'Simpson Index','double',[dimid_lon1 dimid_lat1 dimid_mon1]);
netcdf.putAtt(ncid1,varid_lon1,'standard_name','Longitude');
netcdf.putAtt(ncid1,varid_lon1,'units','Degrees East');
netcdf.putAtt(ncid1,varid_lon1,'reference','Grid Center');
netcdf.putAtt(ncid1,varid_lat1,'standard_name','Latitude');
netcdf.putAtt(ncid1,varid_lat1,'units','Degrees North');
netcdf.putAtt(ncid1,varid_lat1,'reference','Grid Center');
netcdf.putAtt(ncid1,varid_mon1,'standard_name','Month');
netcdf.putAtt(ncid1,varid_mon1,'units','Months Since Dec 1994');
    netcdf.putAtt(ncid1,varid_DS,'standard_name','Simpson Index');
    netcdf.putAtt(ncid1,varid_DS,'units','Unitless');
netcdf.endDef(ncid1)
netcdf.putVar(ncid1,varid_lon1,lon);
netcdf.putVar(ncid1,varid_lat1,lat);
netcdf.putVar(ncid1,varid_mon1,month); 
    netcdf.putVar(ncid1,varid_DS,Simpson);
netcdf.close(ncid1)
clear ncid* dimid* varid*

toc