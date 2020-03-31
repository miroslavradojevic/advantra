clear all; close all; clc;

% run matlab from terminal:
% /usr/local/MATLAB/R2017a/bin/matlab -nodesktop -nosplash -r "foo"

% VLFEAT library
VLFEAT_DIR='/home/miroslav/vlfeat-0.9.20'; % path to the library
run(fullfile(VLFEAT_DIR,'toolbox/vl_setup')); % link vlfeat library
% vl_version verbose % check if library is loaded

ANNOTS_DIR = '/home/miroslav/gndtth/'; % path to the directory with the annotations

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

filterDims = [4 8 16];
%cellSizes = D/[4 8 16]; % [128 64 32];  % 4 8 16
filterSizes = 31 * filterDims.^2;
% 112,142 -- 4x4x31
% 59,66 -- 8x8x31
% 31,32 -- 16x16x31
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tStart = tic;

for lMax = 0 : 1, % scales
for cI = 1 : length(filterDims),
    
    hogFeats = zeros(length(x), 1); % where the computed hog features will be stored (lMax+1)* filterSizes(cI)
    
    for i = 1 : min(Inf, length(x)),

        if i==1 || mosaicIndex(i)~=mosaicIndex(i-1),      
            mosaicName = [ANNOTS_DIR sprintf('m%02d.tif', mosaicIndex(i))];
        end
    
        fprintf(1,'x=%d, y=%d, %s ... ', x(i), y(i), mosaicName);
    
        patch=imread(mosaicName, 'PixelRegion',{[y(i),y(i)+500],[x(i),x(i)+500]});
    
        for l = 0 : lMax,
            
            if l > 0,
                patch1 = imresize(patch, (0.5)^l); % rescale
            else 
                patch1 = patch;
            end
            
            patch1=single(im2double(patch1)*255); % scale to [0,255] and convert to single
            patchSize = size(patch1,1);
            hogImg = vl_hog(patch1, patchSize/filterDims(cI)); % , 'verbose'   cellSizes(cI)
            K = length(hogImg(:));
            
            %fprintf(1, '%d: %d - %d [%d == %d] \n', i, (l*K+1), ((l+1)*K), K, filterSizes(cI));
            
            
            hogFeats(i, (l*K+1):((l+1)*K)) = hogImg(:)';
        
        end
    end
    
    % save the hog features
    fname = strcat('hog', '_', num2str(filterDims(cI)), '_', num2str(lMax));
    
    fprintf(1, 'saving %s... ', fname);
    
    save(strcat(fname, '.mat'), 'hogFeats');
    fprintf(1, 'mat. \n');
    %dlmwrite(strcat(fname, '.csv'), hogFeats, 'delimiter',',', 'precision','%.2f');%too slow to save and too large encoding
    %fprintf(1, 'csv. \n');
       
    clear patch patch1 hogFeats hogImg;
end
end

tElapsed = toc(tStart);
fprintf(1, ' done, elapsed=%.4f hours \n', tElapsed/3600.0);

%[~,d] = vl_sift(patch); % extract sift
