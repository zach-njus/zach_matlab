function [ D ] = span_gap( points,model,curr_idx,bin,labeled,ave_w )
%SPAN_GAP takes in the points that have been previously traversed

pointsc = points(:,1);
pointsr = points(:,2);

%find endpoints that qualify to be jumped to
%{
[rowb,colb] = find(bwmorph(skel,'endpoints'));
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
[distance,~] = min(dist); 
%}

%find endpoints that qualify to be jumped to
[rowb,colb] = find(bwmorph(labeled,'endpoints'));
indexes = [];

for j = 1:length(rowb)
    if(rowb(j) ~= pointsr(end) || colb(j) ~= pointsc(end))
        indecies = pixelLine_pixels([pointsc(end),pointsr(end)],[colb(j),rowb(j)]);
        
        values = zeros(length(indecies),1);
        for k = 1:size(indecies,1)
            values(k) = bin(indecies(k,1),indecies(k,2))
        end
        if(min(values)==1 && length(values)>1)
            indexes = [indexes;j];
        end
    end
end

%if there is a point close enough it can be jumped to 
if(~isempty(indexes))
    %use the old points to predict where the next point should be
    %next_point = predict_point([pointsc(end-15:end),pointsr(end-15:end)],15);
    %use the locaiton of the old model point with the predicted point to
    %find an average location
    %average_point = (next_point+model(curr_idx,:))/2;
    %locs = find(dist<distance*4);
    %locs = find(dist);
    locs = indexes;
    
    %imshow(labeled)
    %hold on;
    %plot(pointsc,pointsr,'r*');
    %plot(colb(locs),rowb(locs),'g*','linewidth',3)
    %plot(average_point(1),average_point(2),'w*','linewidth',2);

    %find the section that appears to match the best
    dist = zeros(length(locs),1);
    dist_sum = dist;
    dist_max = dist;
    num_points = dist;
    for j = 1:length(locs)
        %make temperary model
        temp_model = model;
        %make temperary path
        temp = pixelLine1([colb(locs(j)),rowb(locs(j))],[pointsc(end),pointsr(end)],labeled==labeled(rowb(locs(j)),colb(locs(j))),1);
        temp(pointsr(end),pointsc(end))=1;
        %record current index
        temp_idx = curr_idx;
        %do distance transformation
        D = bwdistgeodesic(temp,pointsc(end),pointsr(end));
        [~,~,~,temp_idx,temp_model] = walk_skel( D,[pointsc,pointsr],temp_model,temp_idx,temp,ave_w );
        dist_idx = zeros(temp_idx-curr_idx,1);
        for k = curr_idx:temp_idx
            [~,dist_idx(k,1)] = vectorRadianDist(temp_model(k,1),temp_model(k,2),model(k,1),model(k,2));
        end
        dist(j,1) = std(dist_idx);
        dist_sum(j,1) = sum(dist_idx);
        dist_max(j,1) = max(dist_idx);
        num_points(j,1) = temp_idx-curr_idx;
        %disp('idxs')
        %disp(temp_idx-curr_idx);
        %disp(dist_sum(j,1)/num_points(j,1))
        %disp(dist_max(j,1));
        %disp(dist(j,1));
    end
    
    %handles small sections having small distances and small standard
    %deviations
    if(sum(num_points<=min(num_points)*2)==1)
        dist(num_points<min(num_points)*2) = 10000;
    end
    [~,loc] = min(dist);
    %close all
    %clc
    %{
    %find the endpoint that is closest to the predicted point
    dist = zeros(length(locs),1);
    for j = 1:length(locs)
        [~,dist(j)] = vectorRadianDist(colb(locs(j)),rowb(locs(j)),average_point(1),average_point(2));
    end
    [~,loc] = min(dist);
    %}
    %connect the current point to the new section
    temp = pixelLine1([colb(locs(loc)),rowb(locs(loc))],[pointsc(end),pointsr(end)],labeled==labeled(rowb(locs(loc)),colb(locs(loc))),1);
    temp(pointsr(end),pointsc(end))=1;
    
    %imshow(temp)
    %hold on;
    %plot(colb(locs),rowb(locs),'g*','linewidth',3)
    %plot(pointsc(end),pointsr(end),'r*');
    %pause(2)
    
    %create a new path to walk along
    D = bwdistgeodesic(temp,pointsc(end),pointsr(end));
    
else
    %there was not a close enough point so the end is reached
    D = nan;
end


end

