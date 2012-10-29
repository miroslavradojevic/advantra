clear all; close all; clc;

current_dir     = pwd;
input_file_dir  = input('folder path? (where ms_convergence_points.m is):  ','s');
cd(input_file_dir);
disp('loading...');
ms_convergence_points;
im = imread('ms_input_image.tif');

% visualize
nr_frames = size(frames, 2);

for iter = 1 : nr_frames,
    
    figure(iter);
    % [row col] is stored
    imshow(im);
    hold on;
    plot(frames{1,iter}(1:2:end,2), frames{1,iter}(1:2:end,1), 'r.');
    export_fig(strcat('frame', num2str(iter),'.tif'), '-tif', '-a1');
    close(iter);
    
end

cd(current_dir);
disp('done.');