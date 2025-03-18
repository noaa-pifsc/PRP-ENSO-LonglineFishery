function [isopleth_depth] = find_isopleth_depth(property_3d_matrix,target_isopleth,depth_max,depth,depth_res_interp)
% Input here are a 3d matrix [depth lat lon] of a property, the target
% isopleth of the property, the maximum depth to be considered, the depth
% vector, and the depth resolution to interpolate onto

[n, o, p] = size(property_3d_matrix); % n is depth levels, o and p are lat and lon
depth_max_loc = dsearchn(depth, depth_max); %We'll limit this exercise to depth_max, and we'll interpolate to every m (depth_res_interp).
isopleth_depth = NaN(o,p);
for q = 1:o
    for r = 1:p
        profile_value = interp1(depth(1:depth_max_loc),squeeze(property_3d_matrix(:,q,r)),1:depth_res_interp:depth_max,'linear'); % linearly interpolate between model depth levels to the levels specified
        profile_depth = 1:depth_res_interp:depth_max;
        if min(profile_value)<target_isopleth  % if the minimm of the profile is less than the target isopleth
            isopleth_depth(q,r) = profile_depth(find(profile_value<target_isopleth,1)); % then return the first location where the value is less than the target isopleth
        elseif isnan(min(profile_value))
            isopleth_depth(q,r) = NaN; % If there are no data, return NaN
        else
%             isopleth_depth(q,r) = NaN;
            isopleth_depth(q,r) = depth_max;  % otherwise return the maximum depth
        end
    end %Close loop for longitues
end %Close loop for latitude

end