clc
clear all
close all

%this script analyzes the results from the tracking software
[file,path] = uigetfile('*.xls');
data = xlsread([path,file]);

%%
%calculate the length
blength = zeros(size(data,2)/2,1);
for i = 1:2:size(data,2)-1
    blength((i+1)/2) = sum(sqrt(sum(transpose(diff([data(:,i),data(:,i+1)]).^2))));
end

plot(blength)
hold on;
horizontal_line = ones(length(blength))*mean(blength);
std_upper = horizontal_line+std(blength);
std_lower = horizontal_line-std(blength);
plot(1:length(blength),horizontal_line,'g-');
plot(1:length(blength),std_lower,'r-');
plot(1:length(blength),std_upper,'r-');
xlabel('Time (frames)');
ylabel('Length (pixels)');
%%
%calculate the velocity
