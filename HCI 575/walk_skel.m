function [ points,labeled,c_length,curr_idx,model ] = walk_skel( D,points,model,curr_idx,labeled,ave_w )
%WALK_SKEL is an algorithm utilizing the bwgeodesic transformation to walk
%along a connected sections of pixels

if(~isempty(points))
    pointsc = points(:,1);
    pointsr = points(:,2);
else
    pointsr = [];
    pointsc = [];
end

%remove nans and inf from the distance matrix
D(isnan(D))=0;
D(isinf(D))=-1;

for j = 1:max(max(D));
    %find the next point to connect to
    [rowc,colc] = find(D == j);
    if(~isempty(rowc))
        pointsr = [pointsr;rowc(1)];
        pointsc = [pointsc;colc(1)];
    
    %eliminate if from the skeleton
    labeled(rowc(1),colc(1)) = 0;
    
    %update the current length of the traversed skeleton
    c_length = sum(sqrt(sum(transpose(diff([pointsc,pointsr]).^2))));
    
    %if the skeleton is has been traversed far enough update the pose of
    %the worm
    if(c_length > (curr_idx-1)*ave_w*2.1 && curr_idx < length(model) && length(pointsr)>1)
        [~,dist] = vectorRadianDist(model(curr_idx,1),model(curr_idx,2),pointsc(end),pointsr(end));
        iterations = round(dist/ave_w*10);
        xdiff = (pointsc(end)-model(curr_idx,1))/iterations;
        ydiff = (pointsr(end)-model(curr_idx,2))/iterations;
        
        %imshow(labeled)
        %hold on;
        %this section is just for outputting purposes
        for k = 1:iterations
            for q = curr_idx:length(model)
            model(q,:) = [model(q,1)+xdiff,...
                model(q,2)+ydiff];
            end
            %model = updatePose(model,curr_idx,1,ave_w*2);
            colors = jet(length(model));
        end
        %{
        for z = 1:curr_idx
            circle(model(z,1),model(z,2),ave_w,colors(z,:));
        end
        %pause(.01)
        %clf('reset')
        %}
        
        %perform one last update to ensure that the model point has been
        %set to the skeleton
        model(curr_idx,:) = [pointsc(end),pointsr(end)];
        model = updatePose(model,curr_idx,1,ave_w*2);
        curr_idx = curr_idx + 1;
        
        %colors = jet(length(model));
        %imshow(second_frame);
        %for k = 1:length(model)
            %circle(model(k,1),model(k,2),ave_w,colors(k,:));
        %end
        %pause(.1)
        %clf('reset')
    end
    end
end
points = [pointsc,pointsr];
%figure;

end

