%Small collection of points to start with
clear all
close all
model = [10,10;100,100;190,190;280,280;370,370;460,460];
modelx = transpose(spline(1:6,transpose(model(:,1)),1:.1:6));
modely = transpose(spline(1:6,transpose(model(:,2)),1:.1:6));
model = [modelx,modely];
axis([1 2500 1 2500])
plot(model(:,1),model(:,2),'b');
plot(model(:,1),model(:,2),'bo');
hold on;
pause(2)
radius = sqrt(90^2+90^2)/10;

%Move a point to another location
index = 51;
xdirection = 20;
ydirection = 20;
writerObj = VideoWriter('Drag Movement.avi');
writerObj.FrameRate = 5;
open(writerObj);
hold off
rad = zeros(51,1);
rad(1:26) = log(400:100:2900);
rad(26:end) = log(2900:-100:400);
leftside = zeros(51,2);
rightside = zeros(51,2);

for i = .1:.05:2*pi
    
    model(index,1) = model(index,1) + xdirection;
    model(index,2) = model(index,2) + ydirection;

    model = updatePose(model,index,1,radius);
    model = updatePose(model,index,-1,radius);
    hold on
    colors = jet(size(model,1));
    for j = 1:size(model,1)
        %plot(model(j,1),model(j,2),'*','color',colors(j,:))
        if(j == 1 || j == size(model,1))
            leftside(j,1) = model(j,1);
            leftside(j,2) = model(j,2);
            rightside(j,1) = model(j,1);
            rightside(j,2) = model(j,2);
        elseif(j>1)
            angle = vectorRadianDist(model(j,1),model(j,2),model(j-1,1),model(j-1,2));
            leftside(j,1) = model(j,1)+cos(angle+pi/2)*rad(j)*2;
            leftside(j,2) = model(j,2)+sin(angle+pi/2)*rad(j)*2;
            rightside(j,1) = model(j,1)+cos(angle-pi/2)*rad(j)*2;
            rightside(j,2) = model(j,2)+sin(angle-pi/2)*rad(j)*2;
        end
    end
    plot(model(:,1),model(:,2),'-g');
    plot(rightside(:,1),rightside(:,2),'-r');
    plot(leftside(:,1),leftside(:,2),'-b');
%     plot(model(:,1),model(:,2),'ro')
    %axis([1 2500 1 2500])
    pause(.1)
    frame = getframe(gcf);
    writeVideo(writerObj,frame);
    %close(gcf)
    clf('reset')
    ydirection = 70*cos(i*5+pi)+5;
    xdirection = 20*abs(cos(i*5))+5;
%     ydirection = 70;
%     xdirection = 20;
end
close(writerObj)

