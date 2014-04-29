function [ skel,bin2,dist_trans ] = brugia_skeleton( frame_n,bin1,background,ave_w )
%BRUGIA_SKELETON creates a skeleton with a brugia worm

%find the difference between the background and current frame
sub2 = uint8(abs(double(background) - double(frame_n)));
sub2 = rgb2gray(sub2);
bin2 = sub2 > 10;

%fill in small holes in the worm
bin2=~bwareaopen(~bin2,100);

%label the second frame and only keep the largest blob with an overlap
%from the previous frame
if(isempty(bin1))
    %label and only record the largest object
    [bin2_labeled,n] = bwlabel(bin2);
    msum=0;
    for j = 1:n
        csum = sum(sum(bin2_labeled==j));
        if(csum>msum)
            bin2 = bin2_labeled==j;
            msum = csum;
        end
    end
else
    [bin2_labeled] = bwlabel(bin2);
    bin2_temp = bin2_labeled.*double(bin1);
    n = unique(bin2_temp);
    msum = 0;
    for j = 1:length(n)
        csum = sum(sum(bin2_labeled==n(j)));
        if(csum>msum && n(j)>0)
            bin2 = bin2_labeled==n(j);
            msum = csum;
        end
    end
end

%develop skeleton for the current frame
dist_trans1 = bwdist(~imfill(bin2,'holes'));
dist_trans = bwdist(~bin2);
dist_trans1 = (dist_trans1>ave_w*.8)&(dist_trans1<ave_w*1.5);
dist_trans = bwmorph(bwmorph((dist_trans>ave_w*.8)&(dist_trans<ave_w*1.5),'thin',inf),'dilate',2);
dist_trans = dist_trans|dist_trans1;

%use thin to get it down to single pixel wide line
skel = bwmorph(dist_trans,'thin',inf);

end

