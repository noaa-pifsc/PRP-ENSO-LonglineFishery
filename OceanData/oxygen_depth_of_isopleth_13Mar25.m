% This script attempts to find the depth of an isopleth of oxygen using the
% GLORYS data served by the Copernicus project.  The directory for the
% GLORYS data needs to be loaded as "folder_name."
% 

% An external function
% called find_isopleth_depth is what does the calculations, and that
% function ingests serveral data and variables, namely one time step of the
% oxygen data, the target isopleth, the maximum depth to consider, the
% depth vector, and the resolution of the vertical profile on which to
% interpolate.  Output are the depths of the isopleth.

% Set the director that is being used
% folder_name = 'D:\Historical_datasets\Copernicus_GLOBAL_MULTIYEAR_BGC_001_029'; %Data are from 1993 through 2024
folder_name = 'Z:\historical_reanalyses\Copernicus_GLOBAL_REANALYSIS_BIO_001_029'; %Data are from 1993 through 2024
year_folders = dir(folder_name); %This is a struct with the information about the directories by year

%extract some of the key paramerters from one file, here just using the
%most recent December
filename = [folder_name '\' year_folders(end).name '\' 'mercatorfreebiorys2v4_global_mean_' year_folders(end).name '12.nc'];
lat = nc_varget(filename,'latitude');
lon = nc_varget(filename,'longitude');
depth = nc_varget(filename,'depth');
% you can use ncread in place of nc_varget if you do not have SNCTOOLS (https://mexcdf.sourceforge.net/) installed.
% lat = ncread(filename,'latitude');
% lon = ncread(filename,'longitude');
% depth = ncread(filename,'depth');

% This is where we can set the domain to be considered.  This is 10-40,
% 180-130 region. "loc" variables here are indices in the original matrices
% that are closest to the specified latitude and longitude limits.
lat_loc_min = dsearchn(lat,10);
lat_loc_max = dsearchn(lat,40);
lon_loc_min = dsearchn(lon,-180);  % This doesn't work if you want to go west of the dateline.  We'll have to treat that circumstance later.
lon_loc_max = dsearchn(lon,-130);

% This is a larger spatial area than defined int he option above.
% lat_loc_min = dsearchn(lat,0);
% lat_loc_max = dsearchn(lat,50);
% lon_loc_min = dsearchn(lon,-180);  % This doesn't work if you want to go west of the dateline.  We'll have to treat that circumstance later.
% lon_loc_max = dsearchn(lon,-100);

% Note that the data run from -180 to 180; not very convenient for us.  We
% could download the whole datasent and shift the grid.  That would save
% time, but that would not be computationally efficient.

% Let's try to create an estimate of the depth of a certain isopleth of oxygen at each model location within a subset domain
target_isopleth = 44.661*2; %original units are mmol/m3 or umol/L.  1 ml/l = 10​00/22.391 = 44.661 µmol/l.  So this corresponds to 2 ml/l.
depth_max = 1000;  % Maximum depth we're considering relevant.  In some locations, the isopleth will be deepther than 1000 m, but we are setting that to 1000.
depth_max_loc = dsearchn(depth, depth_max); %We'll limit this exercise to the first 1000 m, and we'll interpolate to every m
depth_res_interp = 1; % resolution of the interpolated depths will be 1m


valid_folders = [];
for i = 1:length(year_folders)
    if str2num(year_folders(i).name)>1
        valid_folders = [valid_folders i];
    end;
end;

clear o2
o2_isopleth_depth = NaN([length(valid_folders)*12 length(lat_loc_min:lat_loc_max) length(lon_loc_min:lon_loc_max)]); % create a matrix of NaN of appropriate dimensions
o2_isopleth_depth_2mlpl = o2_isopleth_depth;
for i = 1:length(valid_folders)
    clear files_list
    files_list = dir([folder_name '\' year_folders(valid_folders(i)).name '\*.nc']); % list the files in that directory.  Should be 12 per year.
    for j = 1:length(files_list)
        o2 = nc_varget([files_list(j).folder '\' files_list(j).name],'o2',[0 0 lat_loc_min-1 lon_loc_min-1],[1 depth_max_loc length(lat_loc_min:lat_loc_max) length(lon_loc_min:lon_loc_max)]);
        % If SNCTOOLS (https://mexcdf.sourceforge.net/) plugin, can use ncread coupled with permute (next two lines) in place of the line above
        % o2 = ncread([files_list(j).folder '\' files_list(j).name],'o2',flip([1 1 lat_loc_min lon_loc_min]),flip([1 depth_max_loc length(lat_loc_min:lat_loc_max) length(lon_loc_min:lon_loc_max)]));
        % o2 = permute(o2,[3 2 1]); % changes the order of the dimensions to [depth lat lon]
        o2_isopleth_depth_2mlpl((i-1)*12+j,:,:) = find_isopleth_depth(squeeze(o2),target_isopleth,depth_max,depth,depth_res_interp); % it's i-1 because the first year is 1 and we want the count to be 0.
        disp([files_list(j).folder '\' files_list(j).name]); % print out where we are in the processing
    end; %Close loop for months
end; %Close loop for years
[m o p] = size(o2_isopleth_depth_2mlpl);

% animation of the propery
% for i = 1:m
%     pcolor(lon(lon_loc_min:lon_loc_max),lat(lat_loc_min:lat_loc_max),squeeze(o2_isopleth_depth_89(i,:,:)));
%     shading flat;
%     hold on;
%     contour(lon(lon_loc_min:lon_loc_max),lat(lat_loc_min:lat_loc_max),squeeze(o2_isopleth_depth_89(i,:,:)),[400 400],'k');
%     set(gca,'clim',[0 1000]);colorbar;
%     pause(0.1);hold off;
% end;

[o2_isopleth_depth_2mlpl_seasonal_clim, o2_isopleth_depth_2mlpl_seasonal_std] = calculate_climatology(o2_isopleth_depth_2mlpl);
o2_isopleth_depth_2mlpl_std_anom = calculate_std_anomaly(o2_isopleth_depth_2mlpl);
o2_isopleth_depth_2mlpl_anom = calculate_anomaly(o2_isopleth_depth_2mlpl);

% Calculete some anomalies, in case that ends up being useful
subset_latitude = lat(lat_loc_min:lat_loc_max);
subset_longitude = lon(lon_loc_min:lon_loc_max);
subset_o2_2mlpl_depth_1993_2024 = o2_isopleth_depth_2mlpl;
subset_o2_2mlpl_depth_1993_2024_anom = o2_isopleth_depth_2mlpl_anom;
subset_o2_2mlpl_depth_1993_2024_std_anom = o2_isopleth_depth_2mlpl_std_anom;

% Ok, and then save these data.  Save all as a .mat file
save('subset_o2_isopleth_data_1000m_MAX_1993_2024.mat','subset_latitude','subset_longitude','subset_o2_2mlpl_depth_1993_2024','subset_o2_2mlpl_depth_1993_2024_anom','subset_o2_2mlpl_depth_1993_2024_std_anom');


% And write a netcdf, at least for the isopleth data themselves
nc_filename = "O2_2mlpl_depth_qrtdeg.nc";
nccreate(nc_filename,"isopleth_depth","Dimensions",{"longitude",p,"latitude",o,"month",m},"Format","classic"); % necessary to reverse the dimensions here.
nccreate(nc_filename,"month", "Dimensions", {"month",m});
nccreate(nc_filename,"latitude", "Dimensions", {"latitude",o});
nccreate(nc_filename,"longitude", "Dimensions", {"longitude",p});
ncwrite(nc_filename,"latitude",lat(lat_loc_min:lat_loc_max));
ncwrite(nc_filename,"longitude",lon(lon_loc_min:lon_loc_max));
ncwrite(nc_filename,"isopleth_depth",permute(o2_isopleth_depth_2mlpl,[3  2 1]));  % necessary to reverse the dimensions here.
ncwrite(nc_filename,"month",1:m);
ncwriteatt(nc_filename,"longitude","standard_name","longitude");
ncwriteatt(nc_filename,"longitude","units","Degrees east");
ncwriteatt(nc_filename,"latitude","standard_name","latitude");
ncwriteatt(nc_filename,"latitude","units","degrees north");
ncwriteatt(nc_filename,"month","standard_name","month");
ncwriteatt(nc_filename,"month","units","months since Dec 1992");
ncwriteatt(nc_filename,"isopleth_depth","standard_name","depth of 2 ml/l disssolved O2");
ncwriteatt(nc_filename,"isopleth_depth","units","meters");
