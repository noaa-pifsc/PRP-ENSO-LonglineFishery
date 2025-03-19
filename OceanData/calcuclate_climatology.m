function [climatology] = caluclate_climatology(property_3d_matrix,target_isopleth,depth_max,depth,depth_res_interp)
% Input here are a 3d matrix [depth lat lon] of a property, the target
% isopleth of the property, the maximum depth to be considered, the depth
% vector, and the depth resolution to interpolate onto


[n, o, p] = size(property_3d_matrix);
depth_max_loc = dsearchn(depth, depth_max); %We'll limit this exercise to the first 1000 m (depth_max), and we'll interpolate to every m (depth_res_interp).
isopleth_depth = NaN(o,p);
for q = 1:o
    for r = 1:p
        profile = interp1(depth(1:depth_max_loc),squeeze(property_3d_matrix(:,q,r)),1:depth_res_interp:depth_max,'linear');
        if min(profile)<target_isopleth
            isopleth_depth(q,r) = find(profile<target_isopleth,1);
        else
            isopleth_depth(q,r) = NaN;
        end
    end %Close loop for longitues
end %Close loop for latitude

end
