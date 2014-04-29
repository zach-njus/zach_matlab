function [ alpha, magnitude ] = vectorRadianDist( endx,endy,startx,starty )
%Returns radian from first point to last point in range from 0 to >2 pi
%along with the magnitude of the distance form first point to last point
dy = endy-starty;
dx = endx-startx;
magnitude=sqrt(dy^2+dx^2);
    
    if (dx == 0)
        if(dy>0)
            alpha=pi/2;
        elseif(dy<0)
            alpha=3*pi/2;
         elseif(dy == 0)
             alpha = 2*pi;
        end
    else
    
    alpha=atan(dy/dx);
  
    if((dx<0)&&(dy<0))
       alpha = alpha + pi; 
    
    elseif((dx>0)&&(dy<0))
       alpha = alpha + 2*pi; 
    
    elseif((dx<0)&&(dy>0))
       alpha = alpha + pi; 
    
    elseif((dy==0)&&(dx<0))
        alpha=pi;
    end
    
    end

end

