%updatePose_body test

clear all
close all

body(1:20)=points();

for i = 1:20
    body(i).midPoint(1:2)=i*30;
    body(i).leftPoint(1:2)=i*30;
    body(i).rightPoint(1:2)=i*30;
end


body(1).diameter = 0;
body(20).diameter = 0;
for i = 2:19
    if(i<5)
        body(i).diameter = 20;
    elseif(i<15)
        body(i).diameter = 30;
    elseif(i<20)
        body(i).diameter = 20;
    end
end

radius = sqrt(30^2+30^2);
axis([1 2500 1 2500])

%Move a point to another location
index = 10;
xdirection = 200;
ydirection = 200;
%writerObj = VideoWriter('Drag Movement.avi');
%writerObj.FrameRate = 5;
%open(writerObj);
newlocation = points;
oldlocaiton = points;

angles = .1:.05:2*pi;
%for i = .1:.05:2*pi
    
    %body(index).midPoint(1) = body(index).midPoint(1) + xdirection;
    %body(index).midPoint(2) = body(index).midPoint(2) + ydirection;
    newlocation.midPoint(1) = body(index).midPoint(1) + xdirection;
    newlocation.midPoint(2) = body(index).midPoint(2) - ydirection;
    oldlocation.midPoint(1) = body(index).midPoint(1);
    oldlocation.midPoint(2) = body(index).midPoint(2);
    
    body = updatePose_body_v2(body,index,newlocation,-1,radius);
    body(index).midPoint = oldlocation.midPoint;
    body = updatePose_body_v2(body,index,newlocation,1,radius);
    
    body = updatePose_body_v2(body,1,body(1),1,radius);
    hold on
    colors = jet(length(body));
    for j = 1:length(body)
        plot(body(j).midPoint(1),body(j).midPoint(2),'*','color',colors(j,:));
        x = [body(j).leftPoint(1),body(j).rightPoint(1)];
        y = [body(j).leftPoint(2),body(j).rightPoint(2)];
        plot(x,y,'-','color',colors(j,:),'linewidth',3);
    end
    axis([1 2500 1 2500])
    pause(.1)
    %frame = getframe(gcf);
    %writeVideo(writerObj,frame);
    %close(gcf)
    if(i < angles(end))
        clf('reset')
    end
    hold off
    ydirection = 70*cos(i*5+pi)+5;
    xdirection = 20*abs(cos(i*5))+5;
%     ydirection = 70;
%     xdirection = 20;
%end
%close(writerObj)

