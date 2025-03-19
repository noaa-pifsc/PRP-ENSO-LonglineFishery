function [isopleth_depth] = find_isopleth_depth(property_3d_matrix, target_isopleth, depth_max, depth, depth_res_interp)
% Finds the depth at which a given isopleth of a property occurs.

[n, o, p] = size(property_3d_matrix); % n is depth levels, o and p are lat and lon
depth_max_loc = dsearchn(depth, depth_max); % Limit depth search

isopleth_depth = NaN(o, p); % Preallocate output

for q = 1:o
    for r = 1:p
        profile_depth = 1:depth_res_interp:depth_max;
        % ***NOTE*** This line below uses a linear interpolation between model levels.  Could consider using something else.
        profile_value = interp1(depth(1:depth_max_loc), squeeze(property_3d_matrix(:,q,r)), profile_depth, 'linear', NaN); % Handle extrapolation

        % Ensure we have valid data before proceeding
        if all(isnan(profile_value))
            isopleth_depth(q,r) = NaN; % No data available
            continue;
        end
        
        % Find the first depth where the interpolated value drops below the isopleth threshold
        idx = find(profile_value < target_isopleth, 1);  **NOTE*** This was written with a subsurface minimum in mind.  Needs mod for other profile shapes.

        if ~isempty(idx)
            isopleth_depth(q,r) = profile_depth(idx);
        else
            isopleth_depth(q,r) = depth_max; % No crossing found, assign max depth
        end
    end
end

end
