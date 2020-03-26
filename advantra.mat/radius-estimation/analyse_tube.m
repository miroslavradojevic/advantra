clear all;
close all;
clc;

% load it
tube;

%
% figure; plot(rad_dist, val, 'b.'); grid on;
% xlabel('radial dist.');
% ylabel('grey-level');
% title('CYLINDER TUBE CONFIGURATION');

% radius values
r = 2.0: 0.2: max(rad_dist);

s = zeros(length(r), length(val));

for i = 1 : length(r),
    s(i, rad_dist<   r(i)) = -1;
    s(i, rad_dist>=  r(i)) = +1;
end

E  = (s*val')'.*(1./(8*r.*r));


figure;

subplot(412);
plot(r, E, 'r', 'LineWidth', 2);
title('E vs. r');
xlabel('radius');
ylabel('snake Energy'); 
grid on;

% derivative
E_diff = [NaN diff(E, 1)];

% shift it
E_diff_1 = [NaN                 E_diff(2:end)];
E_diff_2 = [NaN NaN             E_diff(3:end)];
E_diff_3 = [NaN NaN NaN         E_diff(4:end)];

th = 4;
consec = 15;
E_shifted = zeros(consec, length(E_diff));
for i = 1 : consec,
    E_shifted(i, :) = [ones(1,i)*NaN E_diff(i+1 : end)];
end

E_diff_final =  [...
    E_shifted; ...
    E_diff];

noNeuriteHere = prod(double(E_diff_final<=th));

subplot(413);
plot(r, E_diff, 'g', 'LineWidth', 2);
title('diff(E) vs. r');
xlabel('radius');
ylabel('snake Energy gradient'); 
hold on;
plot(r,  th*ones(1,length(r)), 'r');
plot(r, -th*ones(1,length(r)), 'r');
grid on;

subplot(411);
plot(rad_dist, val, 'b.'); grid on;
xlabel('radial dist.');
ylabel('voxel grey-level');
title('CYLINDER TUBE CONFIGURATION');

subplot(414);
plot(r, noNeuriteHere, 'm', 'LineWidth', 3);
title('NEURITE');
grid on;
