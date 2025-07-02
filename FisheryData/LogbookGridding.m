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
DeepSets = readtable('DeepSets.csv');


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
                Bigeye(r,c,m) = sum(DeepSets.BIGEYE_KEPT(idx), "omitnan") + ...
                                sum(DeepSets.BIGEYE_RELEASED(idx), "omitnan");
                Swordfish(r,c,m) = sum(DeepSets.SWORDFISH_KEPT(idx), "omitnan") + ...
                                   sum(DeepSets.SWORDFISH_RELEASED(idx), "omitnan");
                Mahi(r,c,m) = sum(DeepSets.MAHI_KEPT(idx), "omitnan") + ...
                              sum(DeepSets.MAHI_RELEASED(idx), "omitnan");
                Yellowfin(r,c,m) = sum(DeepSets.YELLOWFIN_KEPT(idx), "omitnan") + ...
                                   sum(DeepSets.YELLOWFIN_RELEASED(idx), "omitnan");
                Pomfret(r,c,m) = sum(DeepSets.POMFRET_KEPT(idx), "omitnan") + ...
                                 sum(DeepSets.POMFRET_RELEASED(idx), "omitnan");
                Vessels(r,c,m) = length(unique(DeepSets.PERMITNUM(idx)));

                clear idx

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
