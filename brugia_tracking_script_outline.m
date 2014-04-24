%first read in the video
[file,path] = uigetfile('*.mp4');
filename = [path,file];
vid = VideoReader(filename);
%%
%display the first frame to perform segmentation
start_frame = 350;
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
    %find binary worm for the previous frame
    sub = background - first_frame;
    sub = rgb2gray(sub);
    bin1 = sub > 10;
    bin1=~bwareaopen(~bin1,50,4);
    
    %fine binary worm for the current frame
    second_frame = read(vid,i);
    sub2 = background - second_frame;
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
    
    
    flow_image(:,:,2) = skel_broke;
    imshow(flow_image)
    set(gca,'position',[0 0 1 1],'units','normalized')
    
    
    changed_points = [];
    for j = 1:length(worm.center)
        %find initial "chunk" of worm section
        chunk = binary_worm(blank,worm.center(j,:),average_width);
        
        %find overlap in previous frame and current frame
        first_sum = sum(sum(bin1&chunk));
        second_sum = sum(sum(bin2&chunk));
        
        %point needs to be moved if there is an overlap but not a total
        %overlap
        if(abs(first_sum-second_sum)>average_width^2/2)

            %find the angle across the "chunk"
            blob = edge(chunk)&bin2;
            blobe = bwmorph(blob,'endpoints');
            [row,col]=find(blobe);
            
            %if there is two points to form an angle calculate new point
            if(~isempty(col)&&length(col)>1)
                
                %find midpoint of the two intersection points
                midpoint = [mean(col),mean(row)];
                
                %find angle through two intersection points
                [theta,dist]=vectorRadianDist(col(1),row(1),col(2),row(2));
                positive_point = [midpoint(1)+dist/1.9*cos(theta+pi/2),...
                    midpoint(2)+dist/1.9*sin(theta+pi/2)];
                negative_point = [midpoint(1)+dist/1.9*cos(theta-pi/2),...
                    midpoint(2)+dist/1.9*sin(theta-pi/2)];
                
                %draw pixel lines to test which direction (-pi/2, pi/2)
                blob_copy = zeros(size(blob,1),size(blob,2));
                blob_copy = pixelLine1(round(positive_point),round(midpoint),blob_copy,1);
                pos_sum = sum(sum(blob_copy&(chunk&bin2)));
                blob_copy(:,:) = 0;
                blob_copy = pixelLine1(round(negative_point),round(midpoint),blob_copy,1);
                neg_sum = sum(sum(blob_copy&(chunk&bin2)));
                
                %find the amount of overlap included in both positive and
                %negative directions
                dist = average_width*2-(neg_sum+pos_sum);
                
                pchunk = binary_worm(blank,[worm.center(j,1)+dist*cos(theta+pi/2),...
                        worm.center(j,2)+dist*sin(theta+pi/2)],average_width);
                nchunk = binary_worm(blank,[worm.center(j,1)+dist*cos(theta-pi/2),...
                        worm.center(j,2)+dist*sin(theta-pi/2)],average_width);
                nsum = sum(sum(dist_trans&nchunk));
                psum = sum(sum(dist_trans&pchunk));
                
                %decide which direction to move the "chunk"
                if(nsum > psum && nsum > 0)
                    %keep track of changed points
                    changed_points = [changed_points;j];
                    %move "chunk"/ center coordinate of that secion of worm
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
            
        elseif(second_sum > 0)
            %keep track of total overlap points
            changed_points = [changed_points;j];
            circle(worm.center(j,1),worm.center(j,2),average_width,[0,0,1]);
        end
    end
    
    plot(worm.center(changed_points,1),worm.center(changed_points,2),'b-','linewidth',3)
    
    pause(.01);
    %first_frame = second_frame;
    

    
    %find the endpoint closest to one of the ends of the changed points
    %this needs work....
    %starting at the largest section
%    dist = zeros(size(rowe,1),2);
%    for j = 1:size(rowe,1)
%        [~,dist(j,1)] = vectorRadianDist(rowe(j),cole(j),...
%            worm.center(changed_points(1),1),worm.center(changed_points(1),2));
%        [~,dist(j,2)] = vectorRadianDist(rowe(j),cole(j),...
%            worm.center(changed_points(1),1),worm.center(changed_points(1),2));
%    end

    %use geodesic distance transformation to track along skeleton
    D = bwdistgeodesic(skel_broke,cole(biggest_section),rowe(biggest_section));
    D(isnan(D))=0;
    D(isinf(D))=-1;
    pointsr = [];
    pointsc = [];
    for j = 1:max(max(D));
        [rowc,colc] = find(D == j);
        pointsr = [pointsr;rowc(1)];
        pointsc = [pointsc;colc(1)];
    end
    
    %once a branchpoint is reached must decide where to go
    dist = zeros(size(changed_points,1),1);
    vect = dist;
    for j = 1:size(dist,1)
        [vect(j,1),dist(j,1)] = vectorRadianDist(worm.center(changed_points(j),1),...
            worm.center(changed_points(j),2),pointsc(end),pointsr(end));
    end
    [val,loc] = min(dist);
    plot(worm.center(changed_points(loc),1),worm.center(changed_points(loc),2),'c*');
    
    dist = zeros(size(rowb,1),1);
    vect = dist;
    for j = 1:size(rowb,1)
        [vect(j,1),dist(j,1)] = vectorRadianDist(colb(j),rowb(j),pointsc(end),pointsr(end));
    end
    dist(dist==0)=10000;
    [val,loc] = min(dist);
    if(val<20)
        if(sum(dist<(val*2))>1)
            disp('use vect');
        else
            pointsr = [pointsr;rowb(loc)];
            pointsc = [pointsc;colb(loc)];
        end
    else
        disp('I messed up or at the end of the worm');
    end
    
    plot(pointsc,pointsr,'r*')
    
%end
