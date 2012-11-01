clear all; close all; clc;

current_dir     = pwd;
input_file_dir  = input('folder path? (where ms_convergence_points.m is):  ','s');
cd(input_file_dir);
disp('loading...');
ms_log;
disp('loading finished.');
fullscreen = get(0,'ScreenSize');
screen_w = fullscreen(3);
screen_h = fullscreen(4);
im = imread('ms_input_image.tif');
image_h = size(im, 1);
image_w = size(im, 2);

magnification = min(round(screen_h/image_h), round(screen_w/image_w));

if(isdir('frames')), 
    rmdir('frames', 's');
end

if(isdir('msvecs')),
    rmdir('msvecs', 's');
end

if(isdir('hists')),
    rmdir('hists', 's');
end

mkdir frames;
mkdir msvecs;
mkdir hists;

im = imread('ms_input_image.tif');
frames = length(new_pos);
curr_clus_nr = zeros(frames);

[X_plot,Y_plot] = meshgrid(1:size(d{1,1},2), 1:size(d{1,1},1));

for fr = 1 : frames,
    
    figure(fr);
    imshow(im); hold on;
    plot(new_pos{1,fr}(:,2), new_pos{1,fr}(:,1),'b.'); 
    export_fig(strcat('./frames/iter', num2str(fr), '.tif'), '-tif', '-a1', '-native');
    close(fr);
    
    d{1,fr} = sqrt(d{1,fr});
    
    figure(fr);
    mesh(X_plot,Y_plot, d{1,fr} );
    export_fig(strcat('./msvecs/mean_shift_vecs', num2str(fr), '.png'), '-png', '-a1', '-native');%, '-transparent'
    close(fr);
    
    figure(fr);
    [N(fr,:), X(fr,:)] = hist(d{1,fr}(:), 100); 
    bar(X(fr,:),N(fr,:));
    export_fig(strcat('./hists/his', num2str(fr), '.png'), '-png', '-a1', '-native');%, '-transparent'
    close(fr);
    
    curr_clus_nr(fr) = sum(N(fr,:)>0);
    
end  

figure;
plot(curr_clus_nr); 
xlabel('iteration');
ylabel('bins'' number');
grid on;

export_fig('stopping_point.png', '-png', '-a1', '-native', '-transparent');

cd(current_dir);
disp('done.');
