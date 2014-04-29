
clc
close all
clear all
%first read in the video
[file,path] = uigetfile('*.mp4');
filename = [path,file];
vid = VideoReader(filename);

%%
%display the first frame to perform segmentation
start_frame = 1;
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
%make output video
newfile = 'video16';
vido = VideoWriter([newfile,'.mp4'],'MPEG-4');
vido.FrameRate = 10;
open(vido);

worm.center = start_worm.center;

%find the initial worm length
worm_length = sum(sqrt(sum(transpose(diff(worm.center).^2))));

%find binary worm for the first frame
[skel,bin1] = brugia_skeleton( first_frame,[],background,average_width );
    
%array to store all points throughout video
worm_coords = [];

for i = start_frame:start_frame+100; 
    %read in current frame
    second_frame = read(vid,i);
    
    %create a skeleton and segment out the worm
    [skel,bin2,fat_skel] = brugia_skeleton( second_frame,bin1,background,average_width );

    %break apart skeleton
    [skel_broke,biggest_section] = break_skel(skel,average_width);
    
    %label sections
    [labeled,n] = bwlabel(skel_broke,8);

    %find endpoints at the ends of each section
    [rowb,colb] = find(bwmorph(skel_broke,'endpoints'));
    [rowe,cole] = find(bwmorph(skel,'endpoints'));

    %find the size of each section of the skeleton
    labeled_sizes = zeros(n,1);
    for j = 1:n
        labeled_sizes(j) = sum(sum(labeled == j));
    end

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
        %remove nans and inf from the distance matrix
        D(isnan(D))=0;
        D(isinf(D))=-1;

        %walk along the skeleton
        %imshow(second_frame);
        %hold on
        [points,labeled,c_length,current_index,worm.center] = walk_skel( D,[pointsc,pointsr],worm.center,...
            current_index,labeled,average_width );
        pointsc = points(:,1);
        pointsr = points(:,2);
        
        %plot(pointsc,pointsr,'r*');

        %A branch point has been reached so a decision is necessary
        D = span_gap([pointsc,pointsr],worm.center,current_index,bin2,labeled,average_width);


        %if there is no point to jump to break out of the loop
        if(isnan(D))
            break;
        end

        %update current length of the worm
        c_length = sum(sqrt(sum(transpose(diff([pointsc,pointsr]).^2))));
     end

     %update the binary image for the next image
     bin1=bin2;

     %find the current discretized skeleton
     worm.center = skel_handles([pointsc,pointsr],num_points);
     if(flag)
         worm.center = transpose(fliplr(transpose(worm.center)));
     end

     %record the current posture of the worm to be written later
     worm_coords = [worm_coords,worm.center];

     %output the image to show the current track and worm
     colors = jet(length(worm.center));
     imshow(second_frame);
    for j = 1:length(worm.center)
        circle(worm.center(j,1),worm.center(j,2),average_width,colors(j,:));
    end
    waitforbuttonpress;
     %write the shown image to a video
     frame = getframe(gcf);
     writeVideo(vido,frame);
     pause(.1)
     %clf('reset')
end

%write the data to an excel file and close the video
xlswrite([newfile,'.xls'],worm_coords);
close(vido);
    
    

