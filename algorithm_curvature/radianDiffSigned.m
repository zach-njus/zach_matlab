function [diff] = radianDiffSigned(angle1,angle2)

diff=abs(angle1-angle2);

if(diff>pi)
	diff=2*pi-diff;
end

if(angle2>angle1)
    diff=-diff;
end

end


