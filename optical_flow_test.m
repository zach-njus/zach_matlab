% This will compare a binary version of the old worm with a binary version
% of the new worm and then color the sections that need to be updated with
% one color and those that don't with another.

%% read in file
filename = 'SCN_single_close_130327_1.mp4';
vid = VideoReader(filename);

%% create background
background = read(vid,1);
background = rgb2gray(background);
image = background;

mask = roipoly(background);
background = roifill(background,mask);
binary = (background - image)>10;
close all;
%% make old skeleton
imshow(read(vid,1));
[old_skel.x,old_skel.y] = ginput;
pause(.1)

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

%% set up left and right points

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
% figures out locations of left and right points
worm = updatePose_body_v2(worm,1,worm(1),1,ave_dist);
%plots the worm
colors = jet(length(worm));
old_binary = false(size(binary,1),size(binary,2));
for j = 1:length(worm)
    plot(worm(j).midPoint(1),worm(j).midPoint(2),'*','color',colors(j,:));
    x = [worm(j).leftPoint(1),worm(j).rightPoint(1)];
    y = [worm(j).leftPoint(2),worm(j).rightPoint(2)];
    plot(x,y,'-','color',colors(j,:));
    if(j<length(worm))
        old_binary = pixelLine1(round(x(1)),round(y(1)),round(worm(j+1).leftPoint(1)),round(worm(j+1).leftPoint(2)),old_binary,1);
        old_binary = pixelLine1(round(x(2)),round(y(2)),round(worm(j+1).rightPoint(1)),round(worm(j+1).rightPoint(2)),old_binary,1);
    end
end
old_binary = imfill(old_binary,'holes');
%% display the intersection of the old binary and the new binary

intersection = zeros(size(image,1),size(image,2),3);
intersection(:,:,1) = binary;
intersection(:,:,2) = old_binary;
imshow(intersection);


