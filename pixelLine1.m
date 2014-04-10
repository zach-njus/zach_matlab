function [ grayimg ] = pixelLine1( x0,y0,x1,y1,grayimg,intensity )
% Implements the optomized Bresenham's line algorithm to draw a pixel wide
% line
% grayimg should be of type uint8
dx=abs(x1-x0);
dy=abs(y1-y0);
if(x0<x1)
    sx=1;
else
    sx=-1;
end
if(y0<y1)
    sy=1;
else
    sy=-1;
end
err=dx-dy;
while(~((x0==x1)&&(y0==y1)))
    grayimg(y0,x0)=intensity;
    e2=2*err;
    if(e2>-dy)
        err=err-dy;
        x0=x0+sx;
    end
    if(e2<dx)
        err=err+dx;
        y0=y0+sy;
    end
end

end

