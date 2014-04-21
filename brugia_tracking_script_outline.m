%first read in the video
[file,path] = uigetfile('*.avi');
filename = [path,file];
vid = VideoReader(filename);
%%
%display the first frame to perform segmentation
start_frame = 900;
first_frame = read(vid,start_frame);
mask = roipoly(first_frame);
background = first_frame;
for i = 1:3
    background(:,:,i) = roifill(first_frame(:,:,i),mask);
end
%%
%have user highlight centerline
imshow(background - first_frame);
pause(2);
[cols,rows] = ginput;
%%
%discretize user points based on distance array
worm.center = skel_handles([cols,rows],60);
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
%use "optical flow" to decide which sections of the worm need to be updated
flow_image = zeros(size(first_frame,1),size(first_frame,2),3);
blank = zeros(size(first_frame,1),size(first_frame,2));
update = blank;
i=start_frame+1;
%for i = 2:vid.NumberOfFrames
    sub = background - first_frame;
    sub = rgb2gray(sub);
    bin1 = sub > 10;
    bin1=~bwareaopen(~bin1,50,4);
    
    second_frame = read(vid,i);
    sub2 = background - second_frame;
    sub2 = rgb2gray(sub2);
    bin2 = sub2 > 10;
    bin2=~bwareaopen(~bin2,50);
    
    flow_image(:,:,1) = bin1;
    flow_image(:,:,2) = bin2;
    dist_trans1 = bwdist(~imfill(bin2,'holes'));
    dist_trans = bwdist(~bin2);
    dist_trans1 = (dist_trans1>average_width*.8)&(dist_trans1<average_width*1.5);
    dist_trans = bwmorph(bwmorph((dist_trans>average_width*.8)&(dist_trans<average_width*1.5),'thin','inf'),'dilate',1);
    dist_trans=dist_trans|dist_trans1;
    flow_image(:,:,3) = dist_trans;
    imshow(flow_image)
    set(gca,'position',[0 0 1 1],'units','normalized')
    changed_points = [];
    for j = 1:length(worm.center)
        chunk = binary_worm(blank,worm.center(j,:),average_width);
        first_sum = sum(sum(bin1&chunk));
        second_sum = sum(sum(bin2&chunk));
        %point needs to be moved
        if(abs(first_sum-second_sum)>average_width^2/2)
            %circle(worm.center(j,1),worm.center(j,2),average_width,[1,1,1]);
            blob = edge(chunk)&bin2;
            blobe = bwmorph(blob,'endpoints');
            [row,col]=find(blobe);
            if(~isempty(col)&&length(col)>1)
                midpoint = [mean(col),mean(row)];
                [theta,dist]=vectorRadianDist(col(1),row(1),col(2),row(2));
                positive_point = [midpoint(1)+dist/1.9*cos(theta+pi/2),...
                    midpoint(2)+dist/1.9*sin(theta+pi/2)];
                negative_point = [midpoint(1)+dist/1.9*cos(theta-pi/2),...
                    midpoint(2)+dist/1.9*sin(theta-pi/2)];
                blob_copy = zeros(size(blob,1),size(blob,2));
                blob_copy = pixelLine1(round(positive_point),round(midpoint),blob_copy,1);
                pos_sum = sum(sum(blob_copy&(chunk&bin2)));
                blob_copy(:,:) = 0;
                blob_copy = pixelLine1(round(negative_point),round(midpoint),blob_copy,1);
                neg_sum = sum(sum(blob_copy&(chunk&bin2)));
                dist = average_width*2-(neg_sum+pos_sum);
                
                pchunk = binary_worm(blank,[worm.center(j,1)+dist*cos(theta+pi/2),...
                        worm.center(j,2)+dist*sin(theta+pi/2)],average_width);
                nchunk = binary_worm(blank,[worm.center(j,1)+dist*cos(theta-pi/2),...
                        worm.center(j,2)+dist*sin(theta-pi/2)],average_width);
                nsum = sum(sum(dist_trans&nchunk));
                psum = sum(sum(dist_trans&pchunk));
                if(nsum > psum && nsum > 0)
                    %move body by set amount and update rest of body
                    changed_points = [changed_points;j];
                    worm.center(j,:) = [worm.center(j,1)+dist*cos(theta-pi/2),...
                        worm.center(j,2)+dist*sin(theta-pi/2)];
                    circle(worm.center(j,1),worm.center(j,2),average_width,[1,0,1]);
                    %worm.center = updatePose(worm.center,j,1,average_width*2);
                elseif(psum > 0)
                    changed_points = [changed_points;j];
                    worm.center(j,:) = [worm.center(j,1)+dist*cos(theta+pi/2),...
                        worm.center(j,2)+dist*sin(theta+pi/2)];
                    circle(worm.center(j,1),worm.center(j,2),average_width,[1,0,1]);
                    %worm.center = updatePose(worm.center,j,1,average_width*2);
                end
            end
            
        else
            changed_points = [changed_points;j];
            circle(worm.center(j,1),worm.center(j,2),average_width,[0,0,1]);
        end
    end
    
    plot(worm.center(changed_points,1),worm.center(changed_points,2),'b-','linewidth',3)
    
    pause(.01);
    %first_frame = second_frame;
    
%end
