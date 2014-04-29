function [ skel_broke,biggest_section ] = break_skel( skel,average_width )
%BREAK_SKEL takes in a skeleton and removes intersection points

%break apart skeleton by eliminating branch points
branchpoints = bwmorph(skel,'branchpoints');
[rowb,colb]=find(branchpoints);
skel_broke = skel;
for j = 1:length(rowb)
    skel_broke(rowb(j)-2:rowb(j)+2,colb(j)-2:colb(j)+2) = 0;
end  

%find global endpoints
endpoints = bwmorph(skel,'endpoints');
[rowe,cole] = find(endpoints);

%find sections containing endpoints and record which one is the largest
[labeled,~] = bwlabel(skel_broke,8);
endpoint_sections = zeros(length(rowe));
skel_endsections = false(size(skel,1),size(skel,2));
psum = 0;
for j = 1:length(rowe)
    if(labeled(rowe(j),cole(j)))
        endpoint_sections(j) = labeled(rowe(j),cole(j));
        csum = sum(sum(labeled == endpoint_sections(j)));
        skel_endsections = skel_endsections | labeled==endpoint_sections(j);
        if(csum>psum)
            psum = csum;
            biggest_section = j;
        end
    end
end

%eliminate small sections of skeleton and add back in end sections if
%eliminated
skel_broke = bwareaopen(skel_broke,5,8);
skel_broke = skel_broke|skel_endsections;

end

