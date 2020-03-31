function vl_sift_mosaic_v2(D, S, OPOS)
% clc;
% clear all;
% close all;
%--------------------------------------
% take annotated mosaics from the TRAIN_DIR (m01.tif, m01.log;...)
% extract square train patches in each mosaic (balanced training):
% (step 1): sample a regular (xy) grid (step=S*D along x and y) and take those
% that were not overlapping with annotations (less than OPOS of the patch area)
% (step 2): random sample (xy) from those that overlap (more than OPOS of the patch area) with annotations, using
% sampling with replacement and using the amount of overlap with the annotation rectangle as importance
% sampling measure
%--------------------------------------
% take number of frames (keypoint locations) per square patch (DxD) and assign sift descriptor with each
% concatenate all sift descriptors of all frames from all patches from all mosacis in a monster array of sift descriptors
% calculate NWORDS centroids using KMeans method. NWORDS centroids will make the dictionary
% generate codeword histograms for each patch
% by counting the number of frames (therefore histogram) assigned to each codeword (centroid, dictionary element)
% by checking the euclidean distance of the frames sift descriptor towards all the centroids
% and picking the closest centroid for each frame - this is effectively BoW approach
%--------------------------------------
% the sample will have (number of neurons annoted * NPOS) positive patches. careful, it's possible that the set won't be balanced.
%--------------------------------------
VLFEAT_DIR='C:/Users/gamata/Documents/MATLAB/vlfeat-0.9.20'; % vlfeat library path (folder in the same directory as train.m)
run(fullfile(VLFEAT_DIR,'toolbox/vl_setup')); % link vlfeat library
vl_version verbose % check if library is ok
%--------------------------------------
addpath('C:\Users\gamata\Documents\MATLAB\nmachinel.test\export_fig'); % library for quality figure export
%--------------------------------------
addpath('C:\Users\gamata\Documents\MATLAB\nmachinel.test\prtools'); % classification library
prwarning off;
prwaitbar off;
%--------------------------------------
%D=500; % rectangular patch size: width/height
%S=0.5; % defines patch sampling grid step relative to patch size D step=S*D
%OPOS=0.5; % overlap to be considered as positive
%.D.S.OPOS_%d_%3.1f_%3.1f
SIGNATURE=sprintf('D.S.OPOS_%d_%3.1f_%3.1f', D, S, OPOS); %.D.S.OPOS.NW.NPOS_500_0.5_0.5;

OVERLAP_PIX=OPOS*D*D;
IMG_EXT='tif'; % mosaic extension
IMG_EXPORT=1; % save midresults

% check if directory and files exist
TRAIN_DIR='C:/Users/gamata/Pictures/mosaics/train'; % name of the subdir with test tifs/logs m01.tif,m01.log,...
if isdir(TRAIN_DIR) == 0, error('TRAIN directory does not exist'); end
%--------------------------------------
filearray = dir([TRAIN_DIR filesep '*.' IMG_EXT]); % get all files in the directory
NumImages = size(filearray,1); % get the number of images
fprintf(1,'%d images found in %s\n', NumImages, TRAIN_DIR);
if NumImages < 0, error('No image in the directory'); end

%G: I remove to initialize to zeros because there are problems when there are less data.
pdsc=[]; %patch descriptors (concatenate all patch descriptors from all train mosaics)
pcnt=0;  % patch count
pframes=[];

tstart=tic;

datamat=zeros(1,NumImages); %create an array to create NumImages

for i=1:NumImages,
    
    imgname = [TRAIN_DIR filesep filearray(i).name]; % read mosaic
    img=imread(imgname);
    
    %to extract the number of mosaic
    imgName=filearray(i).name;
    [nimg,matches] = strsplit(imgName,{'m','.'},'CollapseDelimiters',true);% only extract the number
    datamat(i)=str2double(nimg(2));%this array saves the numbers of mosaic to be able create the mat files later.
    
    fprintf(1,'%s, [%d, %d] vl_sift():\n', imgname, size(img,1), size(img,2));
    
    % go through the grid, take ones with tagmap=0
    Xgrid=1 : round(S*D) : size(img,1)-D+1;
    Ygrid=1 : round(S*D) : size(img,2)-D+1;
    ppi_all=length(Xgrid)*length(Ygrid); % nr patches per image
    ppi_cnt=0; % patch per image count (to plot progress)
    
    x_patch=zeros(1,length(img(:))); % for vizualizations
    y_patch=zeros(1,length(img(:)));
    c_patch=0;
    
    %to save the coordinates of the patches
    coord=[];%it will have 3 columns: nº mosaic, Xcoord, Ycoord
        
    progress=0;
    tic; % measure time
    
    pframTot=[];%zeros(4,length(img(:)));
    pdscTot=[];
    
    tagmap=zeros(size(img));
    
     intersF=[];
    for x = Xgrid,
        for y = Ygrid,
            % The vl_sift command requires a single precision gray scale image.
            % It also expects the range to be normalized in the [0,255]
            % that's not the case now, patches aren't min-max normalized
            patch=imread(imgname,'PixelRegion',{[x,x+D-1],[y,y+D-1]});
            
            ppi_cnt = ppi_cnt+1;
            progress1=floor(((ppi_cnt/ppi_all)/0.1));
            if progress1~=progress,
                progress=progress1;
                fprintf(1,'%d%% ', progress*10);
            end
         
                c_patch=c_patch+1;
                x_patch(1,c_patch)=x;
                y_patch(1,c_patch)=y;
                
                
                coord=[coord;str2double(nimg(2)) x y];
                %to save the coordenates of the patches to obtain their
                %CHRAM-features later
                auxArray=[str2double(nimg(2)) x y];
                dlmwrite('coordPatches.csv',auxArray,'-append','delimiter',',');
                
                
                % extract sift
                patch=single(im2double(patch)*255); % scale to [0,255] and convert to single
                % max val is not 255 but the max encountered within the patch
                [F,d] = vl_sift(patch); % F: 4 x #frames, d: 128 x #frames
                % add it to the patch descriptor list
                pcnt=pcnt+1;
                pdsc{1,pcnt}=d; % concatenate #frames x 128 matrix
                pframes{1,pcnt}=F;
                auxF=[round(F(1,:)+x); round(F(2,:)+y); F(3,:); F(4,:)];
                auxF=auxF.';
                d=d.';
               
                storeF=[];
                storeD=[];
                for idx=1:length(auxF)
                    
                    if tagmap(auxF(idx,1),auxF(idx,2))==0,
                        tagmap(auxF(idx,1),auxF(idx,2))=1;  
                        storeF=[storeF; auxF(idx,:)];
                        storeD=[storeD; d(idx,:)];
                    else
                        intersF=[intersF; auxF(idx,:)];
                    end
                end
               
               pdscTot=[pdscTot;storeD];
%                auxF=[F(1,:)+x; F(2,:)+y; F(3,:); F(4,:)];
%                auxF=auxF.';
               pframTot=[pframTot;storeF];
            
            
        end % Ygrid
    end % Xgrid
    
     %save the variables per each mosaic in a .mat-file
    auxname=['m0' num2str(datamat(i)) 'mosaic'];
    save([auxname '.mat'], 'coord', 'pdsc','pframes', 'pframTot', 'pdscTot','D');
    
    x_patch=x_patch(1,1:c_patch);
    y_patch=y_patch(1,1:c_patch);
    
     figure, imshow(tagmap, []);
        export_fig(sprintf('%s.tagmap.%s.pdf',filearray(i).name, SIGNATURE), '-pdf', '-a2');
        close all;
    
    fprintf(1, '\n');
    
    %to draw a map (size mosaic) with the keypoints of sift-desccriptors
%     yellow = uint8([255 255 0]); % [R G B]; class of yellow must match class of I
%     shapeInserter = vision.ShapeInserter('Shape','Circles','BorderColor','Custom','CustomBorderColor',yellow);
%     featmap=zeros(size(img));
% %     imshow(featmap);
%     circles = int32([pframTot(:,1) pframTot(:,2) pframTot(:,3)]); %  [x1 y1 radius1;x2 y2 radius2]
%     RGB = repmat(featmap,[1,1,3]); % convert I to an RGB image
%     J = step(shapeInserter, RGB, circles);
% %     figure, imshow(J); 
%     
%     
%    
%     export_fig(sprintf('%s.J.%s.%s.pdf',filearray(i).name,'J',SIGNATURE), '-pdf', '-a2');
%     close all;
   
%     %save the variables per each mosaic in a .mat-file
%     auxname=['m0' num2str(datamat(i)) 'mosaic'];
%     save([auxname '.mat'], 'coord', 'pdsc','pframes', 'pframTot', 'pdscTot','D');
    
    fprintf(1,' %6.2f sec\n',toc);
    
 end
tend=toc(tstart);
fprintf(1,'\n\t%6.2f sec\n',tend);
%-------------------------------------
function ovlp = overlap(x1,y1,h1,w1, x2,y2,h2,w2)
x11=x1;
y11=y1;
x12=x1+w1;
y12=y1+h1;

x21=x2;
y21=y2;
x22=x2+w2;
y22=y2+h2;

xdiff=x22-x11;
if x11<x21,
    xdiff=x12-x21;
end

if xdiff<0,
    xdiff=0;
end

ydiff=y22-y11;
if y11<y21,
    ydiff=y12-y21;
end

if ydiff<0,
    ydiff=0;
end

ovlp=xdiff*ydiff;