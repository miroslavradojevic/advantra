%function spherical_coords
%test the possibilities of projecting small spheres (low resolution) to
% azimuth & polar angle
close all; clear all; clc;
delete layer*.tif
%% define the sphere around (0,0,0)
radius = 10;

range = -2*radius : 2*radius;

x = -2*radius : 2*radius;
y = -2*radius : 2*radius;
z = -2*radius : 2*radius;

len = length(range);

img = zeros(len, len, len);
r   = zeros(len, len, len);


sig = radius/6;

count = 1;

for i = 1 : len,
    for j = 1 : len,
        for k = 1 : len,
            
            %r(i, j, k) = sqrt(x(i)^2+y(j)^2+z(k)^2);
            [ theta, phi, r(i, j, k) ] = cartesian2spherical( x(i), y(j), z(k));

            if r(i, j, k)<=radius,
                
                img(i, j, k) = round(255 * exp(-(r(i, j, k)^2)/(2*sig^2)));
                inclination(count) = theta;
                azimuth(count) = theta;
                
            end
            
            count = count + 1;
            
        end
    end
end

for layer = 1 : len,
    imwrite(img(:,:, layer), sprintf('layer%04d.tif', layer),'tif');
end

% show angles for each of the points
figure;
plot(azimuth, theta, 'kd'); grid on;
disp('done...');




