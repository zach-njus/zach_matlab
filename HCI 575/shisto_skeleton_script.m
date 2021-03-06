
clc
close all
clear all
%first read in the video
[file,path] = uigetfile('*.mp4');
filename = [path,file];
vid = VideoReader(filename);
%%
%display the first frame to perform segmentation
start_frame = 100;
first_frame = read(vid,start_frame);
mask = roipoly(first_frame);
background = first_frame;
for i = 1:3
    background(:,:,i) = roifill(first_frame(:,:,i),mask);
end
%%
%use "optical flow" to decide which sections of the worm need to be updated
flow_image = zeros(size(first_frame,1),size(first_frame,2),3);
average_width=50;
for i = start_frame:start_frame+30;
    %find binary worm for the previous frame
    sub = uint8(abs(double(background) - double(first_frame)));
    sub = rgb2gray(sub);
    bin1 = sub > 10;
    bin1=~bwareaopen(~bin1,50,4);
    
    %fine binary worm for the current frame
    second_frame = read(vid,i);
    sub2 = uint8(abs(double(background) - double(second_frame)));
    flow_image = double(second_frame)/255;
    sub2 = rgb2gray(sub2);
    bin2 = sub2 > 10;
    bin2=~bwareaopen(~bin2,50);
    
    %develop skeleton for the current frame
    dist_trans1 = bwdist(~imfill(bin2,'holes'));
    dist_trans = bwdist(~bin2);
    dist_trans1 = (dist_trans1>average_width*.8)&(dist_trans1<average_width*1.5);
    dist_trans = bwmorph(bwmorph((dist_trans>average_width*.8)&(dist_trans<average_width*1.5),'thin','inf'),'dilate',1);
    dist_trans = dist_trans|dist_trans1;
    
    skel = bwmorph(dist_trans,'thin',inf);
    
    %break apart skeleton
    branchpoints = bwmorph(skel,'branchpoints');
    [rowb,colb]=find(branchpoints);
    skel_broke = skel;
    for j = 1:length(rowb)
        skel_broke(rowb(j)-1:rowb(j)+1,colb(j)-1:colb(j)+1) = 0;
    end  
    branchpoints = bwmorph(skel_broke,'endpoints');
    [rowb,colb]=find(branchpoints);
    
    %find global endpoints
    endpoints = bwmorph(skel,'endpoints');
    [rowe,cole] = find(endpoints);
    
    %find sections containing endpoints
    [labeled,n] = bwlabel(skel_broke,8);
    endpoint_sections = zeros(length(rowe));
    skel_endsections = false(size(skel,1),size(skel,2));
    psum = 0;
    for j = 1:length(rowe)
        if(labeled(rowe(j),cole(j)))
            endpoint_sections(j) = labeled(rowe(j),cole(j));
            csum = sum(sum(labeled == endpoint_sections(j)));
            skel_endsections = skel_endsections | labeled==endpoint_sections(j);
            if(csum>psum)
                psum = csum;
                biggest_section = j;
            end
        end
    end
    
    %eliminate small sections of skeleton and add back in end sections if
    %eliminated
    skel_broke = bwareaopen(skel_broke,10,8);
    skel_broke = skel_broke|skel_endsections;
    [labeled,n] = bwlabel(skel_broke,8);
    
    %find the size of each section of the skeleton
    labeled_sizes = zeros(n,1);
    for j = 1:n
        labeled_sizes(j) = sum(sum(labeled == j));
    end
    
    %flow_image(:,:,2) = flow_image(:,:,2)+double(skel_broke);
    blob = imfill(bwmorph(edge(sub2,'canny',.1)&mask,'dilate',2),'holes');
    flow_image(:,:,2) = flow_image(:,:,2)+double(bwmorph(bwmorph(blob,'thin',inf),'dilate',2));
    flow_image(:,:,3) = flow_image(:,:,3)+double(blob);
    imshow(flow_image)
    set(gca,'position',[0 0 1 1],'units','normalized')
    
    
    imshow(flow_image);
first_frame = second_frame;

    pause(.1)
end
