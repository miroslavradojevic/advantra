clear all; close all; clc;

% VLFEAT library
VLFEAT_DIR='/Users/miroslav/vlfeat-0.9.20'; % path to the library
run(fullfile(VLFEAT_DIR,'toolbox/vl_setup')); % link vlfeat library
% vl_version verbose % check if library is loaded

ANNOTS_DIR = './annotations/'; % path to the directory with the annotations

% ANNOTS_POS = strcat(ANNOTS_DIR, 'pos.csv');
% fid = fopen(ANNOTS_POS);
% POS = textscan(fid, '%d%d%d%d%d%d', 'delimiter', ',');
% fclose(fid);
% clear POS;

ANNOTS_ALL = strcat(ANNOTS_DIR, 'all.txt');
fid = fopen(ANNOTS_ALL);
DATA = textscan(fid, '%d%d%d%d', 'delimiter', ';');
mosaicIndex = DATA{1};
x = DATA{2};
y = DATA{3};
tag = DATA{4};
fclose(fid);
clear DATA fid;

D = 500;
cellSizes = [32 63 127];  % 4 8 16
filterDims = [4 8 16];
% 112,142 -- 4x4x31
% 59,66 -- 8x8x31
% 31,32 -- 16x16x31
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lMax = 1;

for cellSizeIdx = 1 : length(cellSizes),
    %
    hog = []; % where the computed hog features will be stored
    
    
    
    
    % save the hog features
    fname = strcat('hog', '_', num2str(filterDims(cellSizeIdx)), '_', num2str(lMax), '.mat');
    save(fname, 'hog');
    csvwrite(fname, hog);
    
end



%%%%%

for i = 1 : 5, % length(x),

    if i==1 || mosaicIndex(i)~=mosaicIndex(i-1),      
        mosaicName = [ANNOTS_DIR sprintf('m%02d.tif', mosaicIndex(i))];
    end
    
    fprintf(1,'x=%d, y=%d, %s ... ', x(i), y(i), mosaicName);
    
    patch=imread(mosaicName, 'PixelRegion',{[y(i),y(i)+D],[x(i),x(i)+D]});
    
    for l = 0 : lMax,
        if l > 0,
            patch1 = imresize(patch, (0.5)^l); % rescale
        else 
            patch1 = patch;
        end
        patch1=single(im2double(patch1)*255); % scale to [0,255] and convert to single
        hogRow = vl_hog(patch1, cellSize); % , 'verbose'     
    end
    
    hogRow = vl_hog(patch1, cellSize); % , 'verbose'
    hog(i, :) = hogRow();

    
    % add it to the feature output
    
    
    clear patch1 hog;
    
end


%     fprintf(1, ' done, elapsed=%.4f sec \n', tElapsed);
%     tStart = tic;
%     tElapsed = toc(tStart);

    %%%%% sift %%%%%
    %[~,d] = vl_sift(patch); % extract sift
    
    %%%%% hog %%%%%
