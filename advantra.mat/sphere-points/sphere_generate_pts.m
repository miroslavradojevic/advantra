clear all;
close all;
clc;

N = input('number of points? ');
R = input('radius? ');

for k = 1 : N,
    
    h(k)    = -1 + 2 * (k-1)/(N-1); 
    % h = [-1, 1] full sphere 
    % h = [-1, 0], h = [0, 1] half sphere (2 different halfs)
    % range of h determines the sphere size
    theta(k)= acos(h(k));
    if (k==1) || (k==N),
        fi(k)   = 0;
    else
        fi(k)   = mod(fi(k-1)+(3.6/sqrt(N))*(1/sqrt(1-h(k)^2)), 2*pi);
    end
    
    x(k) = R * sin(theta(k)) * cos(fi(k));
    y(k) = R * sin(theta(k)) * sin(fi(k));
    z(k) = R * cos(theta(k));
   
    disp(strcat(num2str(k),' / ', num2str(N)));
    
end

figure;
plot3(x', y', z', 'ro');
grid on;
axis equal;

hold on;
plot3(0, 0, 0, 'bd');

figure;
plot3(x', y', z', 'ro');
grid on;
axis equal;

hold on;
plot3(0, 0, 0, 'bd');

disp('done...');