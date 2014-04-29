classdef points
    %UNTITLED6 Summary of this class goes here
    %   Detailed explanation goes here
    %points are a 1x2 array [x y]
    
    properties
        midPoint;
        leftPoint;
        rightPoint;
        diameter;
    end
    
    methods
        function obj = points(mid,left,right)
            if nargin ~= 0
                obj.midPoint = mid;
                obj.leftPoint = left;
                obj.rightPoint = right;
            else
                obj.midPoint = [0,0];
                obj.leftPoint = obj.midPoint;
                obj.rightPoint = obj.leftPoint;
            end
        end
    end
    
end
%{
classdef DocArrayExample
   properties
      Value
   end
   methods
      function obj = DocArrayExample(F)
         if nargin ~= 0 % Allow nargin == 0 syntax
            m = size(F,1);
            n = size(F,2);
            obj(m,n) = DocArrayExample; % Preallocate object array
            for i = 1:m
               for j = 1:n
                  obj(i,j).Value = F(i,j);
               end
            end
         end
      end
   end
end
%}

