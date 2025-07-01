% The purpose of this script is to make grids of catchability.  
% Catchability here refers to the degree of overlap between species' 
% habitat and longline gear.  It is calculated for each grid cell, with 
% habitat defined was waters with O2 >= 2ml/l and temperatures between 8 
% and 14 degrees Celsius.
%
% Catchability (C) ranges from 0 - 1 and is calculated for each grid cell 
% using equations (1) and (2):
% C = (D_gh)^2 / D_g * D_h			(1)
% D_gh = a - (a - b) - d 		    (2)
% where 
% D_g is maximum gear depth, 
% D_h is the depth range of suitable habitat within a given grid cell, 
% a is the greater of the maximum hook depth or maximum habitat depth, 
% b is the lesser of maximum hook depth or maximum habitat depth, and 
% d is the minimum habitat depth (Woodworth-Jefcoats et al. 2019).  
% Deep-set gear has a maximum gear depth, D_g, of 400 m (Bigelow et al. 
% 2006).  A value of 1 would indicate perfect overlap between a species' 
% vertical habitat and fishing gear and therefore 100% % catchability 
% across the area fished.
%
% For this work, we're using the GODAS and GLORYS data that have been 
% regridded to a common 1-deg x 1-deg grid.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Load temperature data
lat_T = ncread('GODAS_1deg.nc', 'YLAT');
lon_T = ncread('GODAS_1deg.nc', 'XLON');
time_T = ncread('GODAS_1deg.nc', 'TIME');

% Note that the GODAS date units are days since 01-JAN-0001
% whereas Matlab interprets date numbers as days since 00-JAN-0000
% So we're going to add one year (365 days) to GODAS dates
time_vec_T = datevec(time_T + 365); 

depth = ncread('GODAS_1deg.nc', 'LEV1_31');
GODAS = ncread('GODAS_1deg.nc', 'GODAS_REGRID');

% GODAS has the dimensions lon x lat x depth x time
% Rearrange so that lat is dim 1, how Matlab likes things
GODAS = permute(GODAS, [2 1 3 4]); % lat x lon x depth x time

%%% Load oxygen data
lat_O = ncread('O2_2mlpl_depth_1deg_noInterp.nc', 'YLAT');
lon_O = ncread('O2_2mlpl_depth_1deg_noInterp.nc', 'XLON');
time_O = ncread('O2_2mlpl_depth_1deg_noInterp.nc', 'MONTH');
GLORYS = ncread('O2_2mlpl_depth_1deg_noInterp.nc', 'GLORYS_O2_2MLPL_REGRID');

% GLORYS has the dimensions lon x lat x time
% Rearrange so that lat is dim 1, how Matlab likes things
GLORYS = permute(GLORYS, [2 1 3]); % lat x lon x time


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate catchability
% Prepare empty grids that will be used to calculate catchability
a = NaN(length(lat_T), length(lon_T), length(time_T));
b = NaN(length(lat_T), length(lon_T), length(time_T));
d = NaN(length(lat_T), length(lon_T), length(time_T)); 
D_h = NaN(length(lat_T), length(lon_T), length(time_T));

% Define a and b for each grid cell over time - T and O2
for r = 1:1:length(lat_T)
    for c = 1:1:length(lon_T)
        for m = 1:1:length(time_T)
            thermal_hab_idx = find(GODAS(r,c,:,m) >= 8 & GODAS(r,c,:,m) <= 14);
            thermal_hab = depth(thermal_hab_idx);

            % First, check whether there's overlap between thermal habitat
            % and oxygen habit
            % Partial overlap: 2mlpl deeper than 8-14 degrees
            % Select deepest depth with both (limited by temp)
            if max(thermal_hab) <= GLORYS(r,c,m)
                a(r,c,m) = max(400, max(thermal_hab));
                b(r,c,m) = min(400, max(thermal_hab));
                d(r,c,m) = min(thermal_hab);
                D_h(r,c,m) = max(thermal_hab) - min(thermal_hab);
            % Partial overlap: 2mlpl occurs within 8-14 degree range
            % Select deepest depth with both (limited by O2)
            elseif (GLORYS(r,c,m) <= max(thermal_hab) && GLORYS(r,c,m) >= min(thermal_hab))
                a(r,c,m) = max(400, GLORYS(r,c,m));
                b(r,c,m) = min(400, GLORYS(r,c,m));
                d(r,c,m) = min(thermal_hab);
                D_h(r,c,m) = GLORYS(r,c,m) - min(thermal_hab);
            % No overlap
            else
                a(r,c,m) = 0;
                b(r,c,m) = 0;
                d(r,c,m) = 0;
                D_h(r,c,m) = 1; % This is actually zero, 
                                % but it's not mathematically possible
                                % to divide by zero.  Because we set the 
                                % numerator to zero, we can divide by
                                % anything.
            end

        end
    end
end

% Catchability
D_gh = a - (a - b) - d;
C = (D_gh.^2) ./ (400 * D_h);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Package data and save
% Permute the dimensions so that we get (lon, lat, time)
C = permute(C,[2 1 3]);

% Create lat and lon grids
% Horizontally stack column with latitudes
lat_grid = repmat(lat_T,1,size(lon_T,1));
% Transpose column with longitudes to a row and vertically stack
lon_grid = repmat(lon_T',size(lat_T,1),1);

% Save data
% Create file, noting NOT to overwrite an exsiting file
ncid1 = netcdf.create('BigeyeCatchability.nc','NOCLOBBER');

% Define dimensions
dimid_lon1 = netcdf.defDim(ncid1,'Longitude',size(lon_grid,2));
dimid_lat1 = netcdf.defDim(ncid1,'Latitude',size(lat_grid,1));
dimid_mon1 = netcdf.defDim(ncid1,'Month',size(C,3));

% Define variables
varid_lon1 = netcdf.defVar(ncid1,'Longitude','double',dimid_lon1);
varid_lat1 = netcdf.defVar(ncid1,'Latitude','double',dimid_lat1);
varid_mon1 = netcdf.defVar(ncid1,'Month','double',dimid_mon1);

varid_C = netcdf.defVar(ncid1,'Catchability','double',[dimid_lon1 dimid_lat1 dimid_mon1]);

% Define attributes
netcdf.putAtt(ncid1,varid_lon1,'standard_name','Longitude');
netcdf.putAtt(ncid1,varid_lon1,'units','Degrees East');
netcdf.putAtt(ncid1,varid_lon1,'reference','Grid Center');

netcdf.putAtt(ncid1,varid_lat1,'standard_name','Latitude');
netcdf.putAtt(ncid1,varid_lat1,'units','Degrees North');
netcdf.putAtt(ncid1,varid_lat1,'reference','Grid Center');

% netcdf.putAtt(ncid1,varid_z1,'standard_name','Depth');
% netcdf.putAtt(ncid1,varid_z1,'units','meters');
% netcdf.putAtt(ncid1,varid_z1,'note','Data do not have a depth dimension')

netcdf.putAtt(ncid1,varid_mon1,'standard_name','Month');
netcdf.putAtt(ncid1,varid_mon1,'units','Months Since Dec 1994');

netcdf.putAtt(ncid1,varid_C,'standard_name','Catchability');
netcdf.putAtt(ncid1,varid_C,'units','unitless');

netcdf.endDef(ncid1)
% netcdf.reDef(ncid1) in case it's necessary to reenter define mode.

% Put the data in the file
netcdf.putVar(ncid1,varid_lon1,lon_T);
netcdf.putVar(ncid1,varid_lat1,lat_T);
% netcdf.putVar(ncid1,varid_z1,0);
% netcdf.putVar(ncid1,varid_mon1,1:1:size(C,4));
netcdf.putVar(ncid1,varid_mon1,(time_O-24)); % time since Dec 1994 not 1992
netcdf.putVar(ncid1,varid_C,C);

% Close the files so they can be used
netcdf.close(ncid1)
