function [ point ] = predict_point( array, time )
%predict_point takes in an array of points and then calculates the
%acceleration and velocity of the array then makes a prediciton on where a
%future point will be. Array should be organized [X;Y]

velocity = mean(diff(array));
acceleration = mean(diff(diff(array)));

point = zeros(1,2);
point(1,1) = array(end,1)+velocity(1)*time+.5*acceleration(1)*time^2;
point(1,2) = array(end,2)+velocity(2)*time+.5*acceleration(2)*time^2;
end

