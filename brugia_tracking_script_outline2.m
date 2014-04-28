
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

%have user highlight centerline
imshow(uint8(abs(double(background) - double(first_frame))));
pause(2);
[cols,rows] = ginput;
%}
%%
%discretize user points based on distance array

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
newfile = 'video14';
vido = VideoWriter([newfile,'.mp4'],'MPEG-4');
vido.FrameRate = 10;
open(vido);
worm.center = start_worm.center;
%use "optical flow" to decide which sections of the worm need to be updated
flow_image = zeros(size(first_frame,1),size(first_frame,2),3);
i=start_frame+1;
worm_length = sum(sqrt(sum(transpose(diff(worm.center).^2))));

%find binary worm for the previous frame
    sub = uint8(abs(double(background) - double(first_frame)));
    sub = rgb2gray(sub);
    bin1 = sub > 5;
    bin1=~bwareaopen(~bin1,50,4);
    %label and only record the largest object
    [bin1_labeled,n] = bwlabel(bin1);
    csum=0;
    msum=0;
    for j = 1:n
        csum = sum(sum(bin1_labeled==j));
        if(csum>msum)
            bin1 = bin1_labeled==j;
            msum = csum;
        end
    end

    %array to store all points throughout video
    worm_coords = [];
for i = start_frame:start_frame+100;
%for i = 1:vid.NumberOfFrames;
 %i=418;   
    %fine binary worm for the current frame
    second_frame = read(vid,i);
    sub2 = uint8(abs(double(background) - double(second_frame)));
    flow_image = double(second_frame)/255;
    sub2 = rgb2gray(sub2);
    bin2 = sub2 > 10;
    bin2=~bwareaopen(~bin2,100);
     %flow_image = zeros(size(bin1,1),size(bin1,2),3);
 %flow_image(:,:,1) = bin1;
 %flow_image(:,:,3) = bin2;
 %imshow(flow_image)
    %label the second frame and only keep the largest blob with an overlap
    %from the previous frame
    [bin2_labeled] = bwlabel(bin2);
    bin2_temp = bin2_labeled.*double(bin1);
    n = unique(bin2_temp);
    csum = 0;
    msum = 0;
    for j = 1:length(n)
        csum = sum(sum(bin2_labeled==n(j)));
        if(csum>msum && n(j)>0)
            bin2 = bin2_labeled==n(j);
            msum = csum;
        end
    end
    
    %develop skeleton for the current frame
    dist_trans1 = bwdist(~imfill(bin2,'holes'));
    dist_trans = bwdist(~bin2);
    dist_trans1 = (dist_trans1>average_width*.8)&(dist_trans1<average_width*1.5);
    dist_trans = bwmorph(bwmorph((dist_trans>average_width*.8)&(dist_trans<average_width*1.5),'thin',inf),'dilate',2);
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
    skel_broke = bwareaopen(skel_broke,round(average_width*4),8);
    %skel_broke = skel_broke|skel_endsections;
    [labeled,n] = bwlabel(skel_broke,8);
    
    %find endpoints at the ends of each section
    [rowb,colb] = find(bwmorph(skel_broke,'endpoints'));
    
    %find the size of each section of the skeleton
    labeled_sizes = zeros(n,1);
    for j = 1:n
        labeled_sizes(j) = sum(sum(labeled == j));
    end
    
    flow_image(:,:,2) = flow_image(:,:,2)+double(skel_broke);
    imshow(rgb2gray(flow_image).*double(bin2))
    
    %direction
    direction = 1;
    
    %use geodesic distance transformation to track along skeleton
    pointsr = [];
    pointsc = [];
    D = bwdistgeodesic(skel_broke,cole(biggest_section),rowe(biggest_section));
    c_length = 0;
    
    %decide which direction from the old skeleton you are starting
    dist = zeros(2,1);
    flag = 0;
    [~,dist(1,1)] = vectorRadianDist(cole(biggest_section),rowe(biggest_section),...
            worm.center(1,1),worm.center(1,2));
    [~,dist(2,1)] = vectorRadianDist(cole(biggest_section),rowe(biggest_section),...
            worm.center(end,1),worm.center(end,2)); 
    if(dist(2)<dist(1))
        worm.center = transpose(fliplr(transpose(worm.center)));
        flag = 1;
    end

 %start tracking along worm   
 current_index = 1;
 while(c_length < worm_length)   
    D(isnan(D))=0;
    D(isinf(D))=-1;
    for j = 1:max(max(D));
        [rowc,colc] = find(D == j);
        pointsr = [pointsr;rowc(1)];
        pointsc = [pointsc;colc(1)];
        labeled(rowc(1),colc(1)) = 0;
        c_length = sum(sqrt(sum(transpose(diff([pointsc,pointsr]).^2))));
        if(c_length > (current_index-1)*average_width*2.1 && current_index < length(worm.center) && length(pointsr)>1)
            [~,dist] = vectorRadianDist(worm.center(current_index,1),worm.center(current_index,2),pointsc(end),pointsr(end));
            iterations = round(dist/average_width*10);
            xdiff = (pointsc(end)-worm.center(current_index,1))/iterations;
            ydiff = (pointsr(end)-worm.center(current_index,2))/iterations;
            imshow(second_frame);
            hold on
            for k = 1:iterations
                for q = current_index:length(worm.center)
                worm.center(q,:) = [worm.center(q,1)+xdiff,...
                    worm.center(q,2)+ydiff];
                end
                %worm.center = updatePose(worm.center,current_index,1,average_width*2);
                colors = jet(length(worm.center));
            end
            for z = 1:length(worm.center)
                circle(worm.center(z,1),worm.center(z,2),average_width,colors(z,:));
            end
            pause(.01)
            clf('reset')
            
            worm.center(current_index,:) = [pointsc(end),pointsr(end)];
            worm.center = updatePose(worm.center,current_index,1,average_width*2);
            current_index = current_index + 1;
            %colors = jet(length(worm.center));
            %imshow(second_frame);
            %for k = 1:length(worm.center)
                %circle(worm.center(k,1),worm.center(k,2),average_width,colors(k,:));
            %end
            %pause(.1)
            %clf('reset')
        end
        current_label = labeled(pointsr(end),pointsc(end));
    end
    [rowb,colb] = find(bwmorph(skel_broke,'endpoints'));
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
    
    %hold on;
    %plot(pointsc,pointsr,'r*');
    %plot(cole,rowe,'go','linewidth',3)
    
    [val,loc] = min(dist);
    close_points = sum(dist<val*2);
    %is there a single point that is close enough
    %if(close_points<2 && val < 10 && max(dist)>0)
        %plot(colb(loc),rowb(loc),'b*');
        %%D = bwdistgeodesic(temp,colb(loc),rowb(loc));
    if(val<50 && val > 0)
        %analyze pixel values connecting close sections
        locs = find(dist<val*2);
        blank = zeros(size(flow_image,1),size(flow_image,2));
        gray_frame = rgb2gray(second_frame);
        pixel_std = zeros(length(locs),1);
        pixel_min = pixel_std;
        for j = 1:length(locs)
            blank = pixelLine1([pointsc(end),pointsr(end)],[colb(locs(j)),rowb(locs(j))],blank,1);
            [rowt,colt] = find(blank);
            if(length(rowt)<20)
                blank = blank | (labeled(rowb(locs(j)),colb(locs(j)))==labeled);
                [rowt,colt] = find(blank);
            end
            pixels = zeros(length(rowt),1);
            for k = 1:length(rowt)
                pixels(k,1) = gray_frame(rowt(k),colt(k));
            end
            pixel_std(j) = std(pixels);
            pixel_min(j) = min(pixels);

            blank(:,:)=0;
            length(pixels)
        end
        [val,loc]=min(pixel_std);
        %title([num2str(i),'     ',num2str(min(pixel_min))]);
        %if(sum(pixel_std<min(pixel_std)*4)>1)
        if(1)
            next_point = predict_point([pointsc(end-15:end),pointsr(end-15:end)],15);
            average_point = (next_point+3*worm.center(current_index,:))/4;
            imshow(second_frame)
                hold on;
            plot(pointsc,pointsr,'r*');
            plot(colb(locs),rowb(locs),'g*','linewidth',3)
            
            plot(average_point(1),average_point(2),'k*','linewidth',2);
            
            %find the endpoint that is closest to the predicted point
            dist = zeros(length(locs),1);
            for j = 1:length(locs)
                [~,dist(j)] = vectorRadianDist(colb(locs(j)),rowb(locs(j)),average_point(1),average_point(2));
            end
            [val,loc] = min(dist);
            
            title('std was very close');
            %pause(3)
        end
        temp = pixelLine1([colb(locs(loc)),rowb(locs(loc))],[pointsc(end),pointsr(end)],labeled==labeled(rowb(locs(loc)),colb(locs(loc))),1);
        temp(pointsr(end),pointsc(end))=1;
        %imshow(temp)
        %hold on;
        %plot(colb(locs),rowb(locs),'g*','linewidth',3)
        %plot(pointsc(end),pointsr(end),'r*');
        %pause(2)
        D = bwdistgeodesic(temp,pointsc(end),pointsr(end));
        
        %D = bwdistgeodesic(skel_broke,colb(locs(loc)),rowb(locs(loc)));
        %plot(colb(locs(loc)),rowb(locs(loc)),'ko');
    else
        break;
    end
        
    c_length = sum(sqrt(sum(transpose(diff([pointsc,pointsr]).^2))));
 end

 bin1=bin2;

 worm.center = skel_handles([pointsc,pointsr],num_points);
 if(flag)
     worm.center = transpose(fliplr(transpose(worm.center)));
 end
 worm_coords = [worm_coords,worm.center];
 colors = jet(length(worm.center));
 imshow(second_frame);
for j = 1:length(worm.center)
    circle(worm.center(j,1),worm.center(j,2),average_width,colors(j,:));
end
 %for j = 1:length(pointsr)
 %    plot(pointsc(j),pointsr(j),'*','color',colors(j,:));
 %end

 frame = getframe(gcf);
 writeVideo(vido,frame);
 pause(.1)
 %clf('reset')
end
xlswrite([newfile,'.xls'],worm_coords);
close(vido);
    
    

