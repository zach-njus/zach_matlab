%alignment script

close all;
clear all;

vid = VideoReader('Greg Report/SCN_multiple_close_130327_1.mp4');
image = read(vid,1);

image = rgb2gray(image);
mask = roipoly(image);
back = roifill(image,mask);
binary = (back-image)>30;

%imshow(binary);

%imshow(image)
pause(1)
[old_skel.x,old_skel.y] = ginput;
pause(.1)
imshow(image);
x=length(old_skel.x);
old_skel.xs = spline(1:x,old_skel.x,1:.01:x);
old_skel.ys = spline(1:x,old_skel.y,1:.01:x);
dist = zeros(length(old_skel.xs),1);
for i = 2:length(old_skel.xs)
    dist(i) = sqrt((old_skel.xs(i-1)-old_skel.xs(i))^2+(old_skel.ys(i-1)-old_skel.ys(i))^2)+dist(i-1);
end

ave_dist = dist(end)/19;
old_skel.xf = zeros(20);
old_skel.yf = zeros(20);
for i = 1:20
    [~,loc] = min(abs(dist-ave_dist*(i-1)));
    old_skel.xf(i) = old_skel.xs(loc);
    old_skel.yf(i) = old_skel.ys(loc);
end
hold on;
plot(old_skel.xf,old_skel.yf,'*');

%%
skel = bwmorph(binary,'thin',inf);
[endPoints.y,endPoints.x] = find(bwmorph(skel,'endpoints'));
imshow(skel)
hold on
plot(endPoints.x,endPoints.y,'*');

worm(1:20) = points;
for i = 1:20
    worm(i).midPoint = [old_skel.xf(i),old_skel.yf(i)];
    worm(i).leftPoint = [old_skel.xf(i),old_skel.yf(i)];
    worm(i).rightPoint = [old_skel.xf(i),old_skel.yf(i)];
end
worm(1).diameter = 0;
worm(20).diameter = 0;
for i = 2:19
    if(i<5)
        worm(i).diameter = 20;
    elseif(i<15)
        worm(i).diameter = 30;
    elseif(i<20)
        worm(i).diameter = 20;
    end
end
newlocation = points;
newlocation.midPoint = [endPoints.x(1),endPoints.y(1)];
worm = updatePose_body_v2(worm,20,newlocation,-1,ave_dist);
worm(20).rightPoint = worm(20).midPoint;
worm(20).leftPoint = worm(20).midPoint;
%worm = updatePose_body_v2(worm,1,worm(1),1,ave_dist);
imshow(image);
hold on;
colors = jet(length(worm));
for j = 1:length(worm)
    plot(worm(j).midPoint(1),worm(j).midPoint(2),'*','color',colors(j,:));
    x = [worm(j).leftPoint(1),worm(j).rightPoint(1)];
    y = [worm(j).leftPoint(2),worm(j).rightPoint(2)];
    plot(x,y,'-','color',colors(j,:));
end
%%
%rotate end point in -pi/2 to pi/2 from original angle
%theta = vectorRadianDist(worm(19).midPoint(1),worm(19).midPoint(2),worm(20).midPoint(1),worm(20).midPoint(2));

%angles = theta-pi/2:.1:theta+pi/2;
imshow(image);
hold on;
for i = length(worm):-1:1
 worm = align_body( binary,worm,i,-1,ave_dist );
 
 
imshow(image);
hold on;
colors = jet(length(worm));
for j = 1:length(worm)
    plot(worm(j).midPoint(1),worm(j).midPoint(2),'*','color',colors(j,:));
    x = [worm(j).leftPoint(1),worm(j).rightPoint(1)];
    y = [worm(j).leftPoint(2),worm(j).rightPoint(2)];
    plot(x,y,'-','color',colors(j,:));
end
pause(2)
end






