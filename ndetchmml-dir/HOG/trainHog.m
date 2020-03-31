function trainHog()

IMG_EXT='jpg'; % mosaic extension

% addpath('/Users/miroslav/ndethcmml/HOG/'); to call this function

% VLFEAT library
VLFEAT_DIR='/Users/miroslav/vlfeat-0.9.20';
run(fullfile(VLFEAT_DIR,'toolbox/vl_setup')); % link vlfeat library
vl_version verbose % check if library is ok

% clc; clear all; close all;

% check if directory and files exist
% name of the subdir with test jpgs
TRAIN_DIR='/Users/miroslav/ndethcmml.tests/NDET.step.Dout.Opos.Oneg_1.00_256_0.5_0.5/NEURONS';
if isdir(TRAIN_DIR) == 0, error('TRAIN directory does not exist'); end

filearray = dir([TRAIN_DIR filesep '*.' IMG_EXT]); % get all files in the directory

NumImages = size(filearray,1); % get the number of images
fprintf(1,'%d images found in %s\n', NumImages, TRAIN_DIR);
if NumImages < 0, error('No image in the directory'); end


for i=1:1,% NumImages
    imgname = [TRAIN_DIR filesep filearray(i).name]; % read train image
    fprintf(1, '%d - %s | %s\n', i, imgname, filearray(i).name);
    %img=imread(imgname);
    %patch=imread(imgname)%,'PixelRegion',{[x,x+D-1],[y,y+D-1]});
    patch=single(im2double( imread(imgname) ) * 255); % scale to [0,255] and convert to single
    
    
    cellSize = 16;
    hog = vl_hog(patch, cellSize, 'verbose') ;
    imhog = vl_hog('render', hog, 'verbose') ;
    
    clf ; imagesc(imhog) ; colormap gray ;
    set(gca,'dataAspectRatio',[1 1 1])
    
end


