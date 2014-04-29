%% select excel file
[file,path] = uigetfile('*.xls');
filename = [path,file];
%% read in excel data
[data,~,~] = xlsread(filename);

temp = [];
for i = 1:2:size(data,1)-1
    temp = [temp,data(i:i+1,:)];
end

index = find(isnan(temp(1,:)),1,'first');
data = temp(1:2,1:index);
%data(1,:) = temp(1:length(temp)/2);
%data(2,:) = temp(length(temp)/2+1:end);

plot(data(1,:),data(2,:))

%% find curvature of the worm's path
thetaDist = zeros(2,length(data)-1);
kvalPoints = thetaDist(1,1:end-1);

% thetaDist has the angles in the first row and the distances in the
% second row
for i=1:length(thetaDist)
    [thetaDist(1,i),thetaDist(2,i)]=vectorRadianDist(data(1,i),data(2,i),data(1,i+1),data(2,i+1));
end


for i=1:length(kvalPoints)
    if(thetaDist(2,i))
        kvalPoints(1,i)=radianDiffSigned(thetaDist(1,i),thetaDist(1,i+1))/thetaDist(2,i);
    else
        %divide by zero prevention
        kvalPoints(1,i)=radianDiffSigned(thetaDist(1,i),thetaDist(1,i+1))/.01;
    end
end
kvalPoints(isnan(kvalPoints)) = 0;


%twenty point averager
kval_ave = kvalPoints(1,1:end-20);
for i = 1:length(kval_ave)
    kval_ave(i) = mean(kvalPoints(1,(10+i)-9:(10+i)+9));
end

plot(abs(kvalPoints),'b-')

hold on;
plot(abs(kval_ave),'r-');
title([num2str(mean(kvalPoints)),'    ',num2str(mean(kval_ave))]);
xlabel('time'); ylabel('curvature');

%% velocity
%{
Dist = zeros(size(row_data,1)-1,size(row_data,2));
for j = 1:size(row_data,1)-1
    for i=1:size(col_data,2)
        [~,Dist(j,i)]=vectorRadianDist(col_data(j,i),row_data(j,i),col_data(j+1,i),row_data(j+1,i));
    end 
end

pcolor(Dist)
xlabel('index'); ylabel('Time');
title('Velocity Worm 3');
colorbar;
%}