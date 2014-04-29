%% select excel file
[file,path] = uigetfile('*.xlsx');
filename = [path,file];
%% read in excel data
[row_data,~,~] = xlsread(filename,'Frame_Cols');
[col_data,~,~] = xlsread(filename,'Frame_Rows');


%% find curvature of the worm in each frame
thetaDist = zeros(2,size(row_data,2)-1);
kvalPoints = zeros(size(row_data,1),size(thetaDist,2)-1);
for j = 1:size(row_data,1)
    % thetaDist has the angles in the first row and the distances in the
    % second row
    thetaDist = zeros(2,size(row_data,2)-1);
    for i=1:size(col_data,2)-1
        [thetaDist(1,i),thetaDist(2,i)]=vectorRadianDist(col_data(j,i),row_data(j,i),col_data(j,i+1),row_data(j,i+1));
    end

    
    for i=1:size(thetaDist,2)-1
        kvalPoints(j,i)=radianDiffSigned(thetaDist(1,i),thetaDist(1,i+1))/thetaDist(2,i);
    end
end

pcolor(kvalPoints)
xlabel('index'); ylabel('Time');
title('Curvature Worm 3');
colorbar;
figure;
%% velocity
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
