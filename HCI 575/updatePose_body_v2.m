function [ model ] = updatePose_body_v2( model,index,newlocation,direction,radius )
%updatePose is a function that given a model, the particular node moved and
%the radius constraint will update the rest of the model

    if(index >= length(model) && direction > 0)
        model(index).leftPoint=newlocation.midPoint;
        model(index).rightPoint=newlocation.midPoint;
        model(index).midPoint = newlocation.midPoint;
        return;
    elseif(index <= 1 && direction < 0)
        model(index).leftPoint=newlocation.midPoint;
        model(index).rightPoint=newlocation.midPoint;
        model(index).midPoint = newlocation.midPoint;
        return;
    else
        if((index <= 1 && direction > 1) || (index >= length(model) && direction < 1))
                model(index).leftPoint=model(index).midPoint;
                model(index).rightPoint=model(index).midPoint;
        end
        %find angle from index location to the next location indicated by
        %direction
        theta = vectorRadianDist(model(index+direction).midPoint(1),...
            model(index+direction).midPoint(2),newlocation.midPoint(1),...
            newlocation.midPoint(2));
        
        prev_theta = vectorRadianDist(model(index+direction).midPoint(1),...
            model(index+direction).midPoint(2),model(index).midPoint(1),...
            model(index).midPoint(2));
        
        weight = .9*log((length(model)/2-abs(length(model)/2-(index+direction))))/log(length(model)/2);
        if(weight == -inf)
            weight = .1;
        elseif(weight == inf)
            weight = .9;
        end
        %if(index > 3 || index < 17)
        %    weight = .95;
        %end
        prev_theta = prev_theta*weight;
        index
        weight

        theta = theta*(1-weight);
        1-weight
        model(index).midPoint = newlocation.midPoint;
        model(index).leftPoint = newlocation.leftPoint;
        model(index).rightPoint = newlocation.rightPoint;
        
        xdirection = radius*cos((theta+prev_theta));
        ydirection = radius*sin((theta+prev_theta));
        
        newlocation.midPoint(1) = model(index).midPoint(1)+xdirection;
        newlocation.midPoint(2) = model(index).midPoint(2)+ydirection;

        xdirection = model(index+direction).diameter*cos((theta+prev_theta)+pi/2);
        ydirection = model(index+direction).diameter*sin((theta+prev_theta)+pi/2);
        
        newlocation.leftPoint(1)=newlocation.midPoint(1)+xdirection;
        newlocation.rightPoint(1)=newlocation.midPoint(1)-xdirection;
        
        newlocation.leftPoint(2)=newlocation.midPoint(2)+ydirection;
        newlocation.rightPoint(2)=newlocation.midPoint(2)-ydirection;
        

colors = jet(length(model));
for j = 1:length(model)
    plot(model(j).midPoint(1),model(j).midPoint(2),'*','color',colors(j,:));
    hold on;
    x = [model(j).leftPoint(1),model(j).rightPoint(1)];
    y = [model(j).leftPoint(2),model(j).rightPoint(2)];
    plot(x,y,'-','color',colors(j,:));
    title(direction)
end
pause(1)
        model = updatePose_body_v2(model,index+direction,newlocation,direction,radius);
    end






