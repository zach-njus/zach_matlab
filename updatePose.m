function [ model ] = updatePose( model,index,direction,radius )
%updatePose is a function that given a model, the particular node moved and
%the radius constraint will update the rest of the model

    if(index >= size(model,1) && direction > 0)
        return;
    elseif(index <= 1 && direction < 0)
        return;
    else
        %find angle from index location to the next location indicated by
        %direction
        theta = vectorRadianDist(model(index+direction,1),...
            model(index+direction,2),model(index,1),model(index,2));
        xdirection = radius*cos(theta);
        ydirection = radius*sin(theta);
        model(index+direction,1) = model(index,1)+xdirection;
        model(index+direction,2) = model(index,2)+ydirection;

        model = updatePose(model,index+direction,direction,radius);
    end


