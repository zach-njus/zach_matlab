function [ handles ] = skel_handles_pixels( array, num_points )
%SKEL_HANDLES takes in an array of points [x,y] and produces another array
%with num_points. Each of these points will be equal distance apart from
%each other. This distance is given by the total distance divided by the
%number of points.

handles = zeros(num_points,2);

%spline the data
arrays = spline(1:length(array),transpose(array),1:.1:length(array));

dist_points = zeros(length(arrays),1);
for i = 2:length(dist_points)
    dist_points(i,1) = sqrt((arrays(1,i-1)-arrays(1,i))^2+(arrays(2,i-1)-...
        arrays(2,i))^2) + dist_points(i-1);
end

arrays=transpose(arrays);
ave_dist = dist_points(end)/(num_points-1);

for i = 1:length(handles)
    sub = abs(dist_points-ave_dist*(i-1));
    [~,loc] = min(sub);
    handles(i,:) = arrays(loc,:);
end

