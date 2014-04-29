function [ mask ] = draw_point_box( binary,point1,point2 )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

mask = false(size(binary,1),size(binary,2));

mask = pixelLine1(floor(point1.leftPoint(1)),floor(point1.leftPoint(2)),floor(point1.rightPoint(1)),floor(point1.rightPoint(2)),mask,1);
mask = pixelLine1(floor(point2.leftPoint(1)),floor(point2.leftPoint(2)),floor(point2.rightPoint(1)),floor(point2.rightPoint(2)),mask,1);
mask = pixelLine1(floor(point1.leftPoint(1)),floor(point1.leftPoint(2)),floor(point2.leftPoint(1)),floor(point2.leftPoint(2)),mask,1);
mask = pixelLine1(floor(point1.rightPoint(1)),floor(point1.rightPoint(2)),floor(point2.rightPoint(1)),floor(point2.rightPoint(2)),mask,1);
mask = bwmorph(mask,'dilate',1);
mask = imfill(mask,'holes');
end

