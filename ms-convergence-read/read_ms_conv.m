clear all; close all; clc;

current_dir     = pwd;
input_file_dir  = input('folder path? (where ms_convergence_points.m is):  ','s');
cd(input_file_dir);
disp('loading...');
ms_log;
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

mkdir frames;

frames = length(new_pos);
for fr = 1:frames,
    figure(fr);
    imshow(imread('ms_input_image.tif')); hold on;
    if fr==1,
        plot(start_pos(:,2),start_pos(:,1),'b.');
    else
        plot(new_pos{1,fr}(:,2), new_pos{1,fr}(:,1),'b.'); 
        
        %frame_diff = (new_pos{1,fr} - new_pos{1,fr-1});
        %ep = max(sqrt(sum(frame_diff.^2,2)));
        
        export_fig(strcat('./frames/iter', num2str(fr), '.tif'), '-tif', '-a1', '-native');
    end
    
    close(fr);
end    

cd(current_dir);
disp('done.');
% % visualize
% nr_frames = size(frames, 2);
% for iter = 1 : nr_frames,
%     figure(iter);
%     % [row col] is stored
%     imshow(im);
%     hold on;
%     plot(frames{1,iter}(1:2:end,2), frames{1,iter}(1:2:end,1), 'r.');
%     export_fig(strcat('frame', num2str(iter),'.tif'), '-tif', '-a1');
%     close(iter);
% end
