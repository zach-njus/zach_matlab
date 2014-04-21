function [ model ] = align_body( binary,model,index,direction,ave_dist )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

%%this is majorly broken

theta = vectorRadianDist(model(index+direction).midPoint(1),model(index+direction).midPoint(2),model(index).midPoint(1),model(index).midPoint(2));

angles = theta-pi/2:.1:theta+pi/2;
area = zeros(1,length(angles));
for i = 1:length(angles)
    %draw a box and then find the amount of overlap
    temp = model;
    temp(index+direction).midPoint(1) = temp(index).midPoint(1) + ave_dist*cos(angles(i));
    temp(index+direction).midPoint(2) = temp(index).midPoint(2) + ave_dist*sin(angles(i));
    
    temp(index+direction).leftPoint(1) = temp(index+direction).midPoint(1) + temp(index+direction).diameter*cos(angles(i)+pi/2);
    temp(index+direction).leftPoint(2) = temp(index+direction).midPoint(2) + temp(index+direction).diameter*sin(angles(i)+pi/2); 
    
    temp(index+direction).rightPoint(1) = temp(index+direction).midPoint(1) + temp(index+direction).diameter*cos(angles(i)-pi/2);
    temp(index+direction).rightPoint(2) = temp(index+direction).midPoint(2) + temp(index+direction).diameter*sin(angles(i)-pi/2);
    
    mask = draw_point_box(binary,temp(index),temp(index+direction));
    area(i) = sum(sum(mask&binary));
end

%find largest area
[~,loc] = max(area);
area
newlocation = points;

newlocation.midPoint(1) = temp(index).midPoint(1) + ave_dist*cos(angles(loc));
newlocation.midPoint(2) = temp(index).midPoint(2) + ave_dist*sin(angles(loc));

model = updatePose_body_v2(model,index+direction,newlocation,direction,ave_dist);
if(direction>0)
    model = updatePose_body_v2(model,1,model(1),direction,ave_dist);
else
    model = updatePose_body_v2(model,length(model),model(end),direction,ave_dist);
end

%{
if(min(model(index).midPoint == model(index).leftPoint))
    angles = -pi/2:.1:pi/2;
    dist = sqrt((model(index).midPoint(1)-model(index+direction).midPoint(1))^2+...
        (model(index).midPoint(2)-model(index).midPoint(2+direction))^2);
    areas = zeros(length(angles));
    for i = 1:length(angles)
        mask = false(size(binary,1),size(binary,2));
        temp_model = model;
        newPoint = zeros(1,2);
        newPoint(1,1) = temp_model(index).midPoint(1)+dist*cos(angles(i));
        newPoint(1,2) = temp_model(index).midPoint(2)+dist*sin(angles(i));
        temp_model(index+direction).midPoint = newPoint;
        
        p1=temp_model(2).leftPoint;
        p2=temp_model(2).rightPoint;
        p3=temp_model(2+direction.leftPoint);
        p4=temp_model(index+direction.rightPoint);
        %mask = pixelLine1(
    end
    
end
%}
end

