%{
clc
close all
clear all
%first read in the video
[file,path] = uigetfile('*.mp4');
filename = [path,file];
vid = VideoReader(filename);
%%
%display the first frame to perform segmentation
start_frame = 300;
first_frame = read(vid,start_frame);
mask = roipoly(first_frame);
background = first_frame;
for i = 1:3
    background(:,:,i) = roifill(first_frame(:,:,i),mask);
end
%%

%have user highlight centerline
imshow(uint8(abs(double(background) - double(first_frame))));
pause(2);
[cols,rows] = ginput;
%}
%%
%discretize user points based on distance array
start_frame = 500;
num_points = 60;
worm.center = skel_handles([cols,rows],num_points);
start_worm.center = worm.center;
imshow(first_frame)
hold on
plot(worm.center(:,1),worm.center(:,2),'r*');
%%
%develop body of the worm
average_width = sqrt((worm.center(1,2)-worm.center(2,2))^2+(worm.center(1,1)-worm.center(2,1))^2)/2;
colors = jet(length(worm.center));
for i = 1:length(worm.center)
    circle(worm.center(i,1),worm.center(i,2),average_width,colors(i,:));
end

%%
vido = VideoWriter('pixel_intensity1.mp4');
vido.FrameRate = 10;
open(vido);
worm.center = start_worm.center;
%use "optical flow" to decide which sections of the worm need to be updated
flow_image = zeros(size(first_frame,1),size(first_frame,2),3);
i=start_frame+1;
worm_length = sum(sqrt(sum(transpose(diff(worm.center).^2))));
for i = start_frame:start_frame+100;
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
    
    %find endpoints at the ends of each section
    [rowb,colb] = find(bwmorph(skel_broke,'endpoints'));
    
    %find the size of each section of the skeleton
    labeled_sizes = zeros(n,1);
    for j = 1:n
        labeled_sizes(j) = sum(sum(labeled == j));
    end
    
    flow_image(:,:,2) = flow_image(:,:,2)+double(skel_broke);
    %figure;
    %imshow(flow_image)
    %set(gca,'position',[0 0 1 1],'units','normalized')
    
    
    pause(.01);
    %figure;
    imshow(rgb2gray(flow_image))
    
    
    %direction
    direction = 1;
    
    %use geodesic distance transformation to track along skeleton
    pointsr = [];
    pointsc = [];
    D = bwdistgeodesic(skel_broke,cole(biggest_section),rowe(biggest_section));
    c_length = 0;
 while(c_length < worm_length)   
    D(isnan(D))=0;
    D(isinf(D))=-1;
    for j = 1:max(max(D));
        [rowc,colc] = find(D == j);
        pointsr = [pointsr;rowc(1)];
        pointsc = [pointsc;colc(1)];
        current_label = labeled(pointsr(end),pointsc(end));
    end
    
    dist = zeros(length(rowb),1);
    vect = dist;
    for j = 1:length(rowb)
        [vect(j),dist(j)] = vectorRadianDist(colb(j),rowb(j),pointsc(end),pointsr(end));
        %prevention from jumping back to previous section
        if(labeled(rowb(j),colb(j))==labeled(pointsr(end),pointsc(end)))
            dist(j) = 0;
        end
    end
    dist(dist==0) = max(dist);
    
    hold on;
    plot(pointsc,pointsr,'r*');
    
    [val,loc] = min(dist);
    close_points = sum(dist<val*2);
    %is there a single point that is close enough
    if(close_points<2 && val < 10)
        plot(colb(loc),rowb(loc),'b*');
        D = bwdistgeodesic(skel_broke,colb(loc),rowb(loc));
    else
        %analyze pixel values connecting close sections
        locs = find(dist<val*2);
        blank = zeros(size(flow_image,1),size(flow_image,2));
        gray_frame = rgb2gray(second_frame);
        pixel_std = zeros(length(locs),1);
        for j = 1:length(locs)
            blank = pixelLine1([pointsc(end),pointsr(end)],[colb(locs(j)),rowb(locs(j))],blank,1);
            [rowt,colt] = find(blank);
            pixels = zeros(length(rowt),1);
            for k = 1:length(rowt)
                pixels(k,1) = gray_frame(rowt(k),colt(k));
            end
            pixel_std(j) = std(pixels);
            %imshow(double(gray_frame)/255+blank);
            %pause(3);
            blank(:,:)=0;
        end
        [val,loc]=min(pixel_std);
        
        D = bwdistgeodesic(skel_broke,colb(locs(loc)),rowb(locs(loc)));
        plot(colb(locs(loc)),rowb(locs(loc)),'ko');
        %{
        %is there a predicted point that is very close
        pre_val = val;
        weight = .1;
        while(weight<1.1)
            point = predict_point([pointsc(end-20:end),pointsr(end-20:end)],pre_val*weight);
            plot(point(1),point(2),'co');

            dist = zeros(length(rowb),1);
            vect = dist;
            for j = 1:length(rowb)
                [vect(j),dist(j)] = vectorRadianDist(colb(j),rowb(j),point(1),point(2));
                %prevention from jumping back to previous section
                if(labeled(rowb(j),colb(j))==labeled(pointsr(end),pointsc(end)))
                    dist(j) = 0;
                end
            end
            dist(dist==0) = max(dist);
            [val,loc] = min(dist);
            close_points = sum(dist<val*1.5);
            %if(close_points > 1)
                weight = weight+.1;
            %else
                %break;
            %end
        end
        
        if(close_points<2 && val < 10)
            D = bwdistgeodesic(skel_broke,colb(loc),rowb(loc));
            plot(colb(loc),rowb(loc),'ko');
        elseif(close_points>2 && val < 10)
            %if there are two points that are the same distance then choose
            %the bigger one 
            locs = find(dist<val*2);
            sizes = zeros(length(locs),1);
            for j = 1:length(locs)
                sizes(j) = sum(sum(labeled==labeled(rowb(locs(j)),colb(locs(j)))));
            end
            locsizes = find(max(sizes)==sizes);
            if(length(locsizes)>1)
                disp('segments are equally spaced and the same size :(');
            else
                D = bwdistgeodesic(skel_broke,colb(locs(locsizes)),rowb(locs(locsizes)));
                plot(colb(locs(locsizes)),rowb(locs(locsizes)),'ko');
            end
        else
            
            break;
        end
        %}
    end
        
    c_length = sum(sqrt(sum(transpose(diff([pointsc,pointsr]).^2))));
 end
 
 colors = jet(length(pointsr));
 for j = 1:length(pointsr)
     plot(pointsc(j),pointsr(j),'*','color',colors(j,:));
 end
 frame = getframe(gcf);
 writeVideo(vido,frame);
 pause(1)
 clf('reset')
end
close(vido);
    
    

