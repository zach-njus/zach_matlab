function [ model ] = updatePose_body( model,index,direction,radius )
%updatePose is a function that given a model, the particular node moved and
%the radius constraint will update the rest of the model

    if(index >= length(model) && direction > 0)
        model(index).leftPoint=model(index).midPoint;
        model(index).rightPoint=model(index).midPoint;
        disp('done')
        return;
    elseif(index <= 1 && direction < 0)
        model(index).leftPoint=model(index).midPoint;
        model(index).rightPoint=model(index).midPoint;
        return;
    else
        if((index <= 1 && direction > 1) || (index >= length(model) && direction < 1))
                model(index).leftPoint=model(index).midPoint;
                model(index).rightPoint=model(index).midPoint;
        end
        %find angle from index location to the next location indicated by
        %direction
        theta = vectorRadianDist(model(index+direction).midPoint(1),...
            model(index+direction).midPoint(2),model(index).midPoint(1),...
            model(index).midPoint(2));
        
        xdirection = radius*cos(theta);
        ydirection = radius*sin(theta);
        
        model(index+direction).midPoint(1) = model(index).midPoint(1)+xdirection;
        model(index+direction).midPoint(2) = model(index).midPoint(2)+ydirection;
        
        xdirection = model(index+direction).diameter*cos(theta+pi/2);
        ydirection = model(index+direction).diameter*sin(theta+pi/2);
        
        model(index+direction).leftPoint(1)=model(index+direction).midPoint(1)+xdirection;
        model(index+direction).rightPoint(1)=model(index+direction).midPoint(1)-xdirection;
        
        model(index+direction).leftPoint(2)=model(index+direction).midPoint(2)+ydirection;
        model(index+direction).rightPoint(2)=model(index+direction).midPoint(2)-ydirection;
        model = updatePose_body(model,index+direction,direction,radius);
    end




