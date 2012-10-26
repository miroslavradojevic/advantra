clear all; close all; clc;

current_dir = pwd;
fold = input('folder path? (where ms_convergence_points.m is):  ','s');
cd(fold);% /home/miroslav/mean_shift_26-10-2012-16_44_44/;
disp('loading...');
ms_convergence_points;
cd(current_dir);

% visualize
nr_frames = size(conv_points{1,1}, 1);

rows = cell(1,nr_frames);
cols = cell(1,nr_frames);

for iter = 1 : nr_frames,
    
    rows{iter} = zeros(numel(conv_points),1);
    cols{iter} = zeros(numel(conv_points),1);
    
    for i = 1 : numel(conv_points),
        rows{iter}(i,1) = conv_points{}(iter,1);
        cols{iter}(i,1) = conv_points{}(iter,2);
    end
end

disp('done.');