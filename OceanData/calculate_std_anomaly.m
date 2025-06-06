function [std_anomaly] = calculate_std_anomaly(property_3d_matrix)
% Input here is a 4d matrix [months lat lon] of a property

[m, o, p] = size(property_3d_matrix);

len = m;
climatology_mean = NaN([12 o p]);
climatology_std = climatology_mean;
for i = 1:12
    climatology_mean(i,:,:) = mean(property_3d_matrix(i:12:len-mod(len,12),:,:));
    climatology_std(i,:,:) = std(property_3d_matrix(i:12:len-mod(len,12),:,:));
end

std_anomaly = NaN([m o p]);
for i = 1:12
    std_anomaly(i:12:end,:,:) = (property_3d_matrix(i:12:end,:,:)-repmat(climatology_mean(i,:,:),[len/12 1 1]))./repmat(climatology_std(i,:,:),[len/12 1 1]);
end

end
