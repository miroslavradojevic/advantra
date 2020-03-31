clear all; close all; clc;

filterDims = 4; % 4 8 16
lMax = 1; % 0 1

% run matlab from terminal:
% /usr/local/MATLAB/R2017a/bin/matlab -nodesktop -nosplash -r "magic(5)"

% VLFEAT library
VLFEAT_DIR='/home/miroslav/vlfeat-0.9.20'; % path to the library
run(fullfile(VLFEAT_DIR,'toolbox/vl_setup')); % link vlfeat library

ANNOTS_DIR = '/home/miroslav/gndtth/'; % path to the directory with the annotations

ANNOTS_ALL = strcat(ANNOTS_DIR, 'all.txt');
fid = fopen(ANNOTS_ALL);
DATA = textscan(fid, '%d%d%d%d', 'delimiter', ';');
mosaicIndex = DATA{1};
x = DATA{2};
y = DATA{3};
tag = DATA{4};
fclose(fid);
clear DATA fid;

% 128 64 32  % 4 8 16
% 112,142 -- 4x4x31
% 59,66 -- 8x8x31
% 31,32 -- 16x16x31
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tStart = tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    hogFeats = zeros(length(x), 1); % where the computed hog features will be stored
    
    for i = 1 : length(x),

        %if i==1 || mosaicIndex(i)~=mosaicIndex(i-1),      
            %
        %end
        
        mosaicName = [ANNOTS_DIR sprintf('m%02d.tif', mosaicIndex(i))];
    
        iRow = y(i)+1;
        iCol = x(i)+1;
        patch=imread(mosaicName, 'PixelRegion',{[iRow, iRow+499],[iCol, iCol+499]});
        
        % pad it with zeros in case it is not [500,500]
        padRow = 500 - size(patch,1);
        padCol = 500 - size(patch,2);
        
        fprintf(1,'%d/%d (%.2f), x=%d-%d, y=%d-%d, %s  patch[%dx%d] pad=[%d %d] \n',... 
        i, length(x), (i*1.0)/length(x), x(i), x(i)+499, y(i), y(i)+499, mosaicName, size(patch,1), size(patch,2),...
        padRow, padCol);
        
            
        if (padRow>0 || padCol>0),
            patch = padarray(patch, [padRow padCol],'replicate','post');
        end
    
        for l = 0 : lMax,
            
            if l > 0,
                patch1 = imresize(patch, (0.5)^l); % rescale
            else 
                patch1 = patch;
            end
            
            patch1=single(im2double(patch1)*255); % scale to [0,255]?! and convert to single
            patchSize = size(patch1,1);
            hogImg = vl_hog(patch1, patchSize/filterDims);
            
            K1 = filterDims^2*31;
            hogFeats(i, (l*K1+1):((l+1)*K1)) = hogImg(:)'; % should fit the range
        
        end
    end
    
    % save the hog features
    fname = strcat('hog', '_', num2str(filterDims), '_', num2str(lMax));
    
    fprintf(1, 'saving %s... ', fname);
    
    save(strcat(fname, '.mat'), 'hogFeats', 'mosaicIndex', 'tag');
    fprintf(1, 'mat. \n');
       
    clear patch patch1 hogFeats hogImg x y tag mosaicIndex;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tElapsed = toc(tStart);
fprintf(1, ' done, elapsed=%.4f hours \n', tElapsed/3600.0);

