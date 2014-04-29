% Just doing background subtraction for each frame to visualize overlapping
% effect (optical flow)

filename = 'shisto.avi';
vid = VideoReader(filename);

background = read(vid,1);
background = rgb2gray(background);
mask = roipoly(background);
background = roifill(background,mask);

intersection = zeros(size(background,1),size(background,2),3);
for i = 2:vid.NumberOfFrames
    curr_image = read(vid,i);
    orig_image = curr_image;
    curr_image = rgb2gray(curr_image);
    intersection(:,:,1) = intersection(:,:,3);
    intersection(:,:,3) = double((background - curr_image)>10);
    %imshow(double(orig_image)/255+intersection);
    imshow(intersection);
    pause(.1)
end