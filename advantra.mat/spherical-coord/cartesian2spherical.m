function [ theta, phi, radius ] = cartesian2spherical( x, y, z )
%CONVERSION TO SPHERICAL COORDINATES 
%   theta - inclination [0, pi]
%   phi   - azimuth     [0, 2*pi)
%   radius 

radius = sqrt(x^2 + y^2 + z^2);

if(radius==0),
    theta = 0;
else
    theta = acos(z/radius);
end

phi = atan2(y, x);

end

