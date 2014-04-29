function [ image ] = binary_worm( image,  center, radius )
%BINARY_WORM this funciton creates a mask of the worm in binary form.
temp = image;
i=1;
%for i = 1:length(center)
    ang=0:0.1:2*pi; 
    xp=radius*cos(ang);
    yp=radius*sin(ang);
    hold on;
    x = center(i,1);
    y = center(i,2);
    for i = 1:length(xp)
        temp(round(y+yp(i)),round(x+xp(i))) = 1;
    end
    image = image | imfill(temp,'holes');
%end

end

