function vl_sift_mosaic(D, S, OPOS, NPOS)
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
%NW=60; % number of clusters, codewords to form the dictionary with
%NPOS=33; % number of patches for each patch marked like neuron by Miguel.
%.D.S.OPOS.NW.NPOS_%d_%3.1f_%3.1f_%d_%d
SIGNATURE=sprintf('D.S.OPOS.NW.NPOS_%d_%3.1f_%3.1f_%d_%d', D, S, OPOS, NW, NPOS); %.D.S.OPOS.NW.NPOS_500_0.5_0.5_75_33;

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

%G: I remove to initialize to zeros because there are problems is there are less data.
pdsc=[]; %patch descriptors (concatenate all patch descriptors from all train mosaics)
ptag=[]; %patch class tag (0=BACKGROUND,1=NEURON)
pcnt=0;  % patch count

tstart=tic;

datamat=zeros(1,NumImages); %create an array to create NumImages

for i=1:NumImages,
    
    imgname = [TRAIN_DIR filesep filearray(i).name]; % read train mosaic
    img=imread(imgname);
    
    %to extract the number of mosaic
    imgName=filearray(i).name;
    [nimg,matches] = strsplit(imgName,{'m','.'},'CollapseDelimiters',true);% only extract the number
    datamat(i)=str2double(nimg(2));%this array saves the numbers of mosaic to be able create the mat files later.
    
    
    annotname = strcat(imgname(1:end-3), 'log'); % read LOG with annotations
    
    neurons=zeros(1000,4); % accumulate rows that correspond to the rectangles [x,y,w,h], allocate 1000 per mosaic
    neuronscnt=0;
    fileID = fopen(annotname);
    while (~feof(fileID)),
        C = textscan(fileID,'%s %d %d %d %d\n');% NEURON COL ROW D D (x=COL, y=ROW in imagej)
        neuronscnt=neuronscnt+1;
        neurons(neuronscnt,:)=[C{2} C{3} C{4} C{5}]; % COL ROW D D
    end
    fclose(fileID);
    neurons=neurons(1:neuronscnt,:);% crop unused allocations
    
    % generate binary tagmap, location map and weight map for the location
    % overlapping with the annotations
    % location map and weight map used for sampling positives in (2)
    tagmap=zeros(size(img)); % to discard positives while going through the grid
    xpos=zeros(1,length(img(:)));
    ypos=zeros(1,length(img(:)));
    wpos=zeros(1,length(img(:)));
    cpos=0;
    for j=1:neuronscnt, % go through the annotated neurons
        x1=neurons(j,2)-neurons(j,4); % ROW
        x2=neurons(j,2)+neurons(j,4);
        y1=neurons(j,1)-neurons(j,3); % COL
        y2=neurons(j,1)+neurons(j,3);
        for x=x1:x2, % go through the neighbourhood of each annot
            for y=y1:y2,
                if x>0 && x<=size(img,1) && y>0 && y<=size(img,2),
                    ovlp = overlap(x,y,D,D, neurons(j,2),neurons(j,1),neurons(j,4),neurons(j,3));
                    if ovlp>=OVERLAP_PIX,
                        tagmap(x,y)=1;
                        cpos = cpos+1;
                        xpos(1,cpos) = x;
                        ypos(1,cpos) = y;
                        wpos(1,cpos) = ovlp/(D*D);
                    end
                end
            end
        end
    end
    
    if cpos==0, error('No positives among mosaic patches'); end
    
    xpos=xpos(1,1:cpos);
    ypos=ypos(1,1:cpos);
    wpos=wpos(1,1:cpos);
    
    fprintf(1,'%s, [%d, %d] vl_sift():\n', imgname, size(img,1), size(img,2));
    
    % go through the grid, take ones with tagmap=0
    Xgrid=1 : round(S*D) : size(img,1)-D+1;
    Ygrid=1 : round(S*D) : size(img,2)-D+1;
    ppi_all=length(Xgrid)*length(Ygrid); % nr patches per image
    ppi_cnt=0; % patch per image count (to plot progress)
    
    xneg=zeros(1,length(img(:))); % for vizualizations
    yneg=zeros(1,length(img(:)));
    cneg=0;
    
    %to save the coordinates of the patches
    coord=[];%it will have 3 columns: nº mosaic, Xcoord, Ycoord
        
    progress=0;
    tic; % measure time
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
            
            if tagmap(x,y)==0,
                % xneg, yneg
                cneg=cneg+1;
                xneg(1,cneg)=x;
                yneg(1,cneg)=y;
                
                
                coord=[coord;str2double(nimg(2)) x y];
                %to save the coordenates of the patches to obtain their
                %CHRAM-features later
                auxArray=[str2double(nimg(2)) x y];
                dlmwrite('coordPatches.csv',auxArray,'-append','delimiter',',');
                
                
                % extract sift
                patch=single(im2double(patch)*255); % scale to [0,255] and convert to single
                % max val is not 255 but the max encountered within the patch
                [~,d] = vl_sift(patch); % f(not used): 4 x #frames, d: 128 x #frames
                % add it to the patch descriptor list
                pcnt=pcnt+1;
                pdsc{1,pcnt}=d; % concatenate #frames x 128 matrix
                ptag(1,pcnt)=0; % class tag
            end
            
        end % Ygrid
    end % Xgrid
    
    xneg=xneg(1,1:cneg);
    yneg=yneg(1,1:cneg);
    
    if cneg==0, error('No negatives among mosaic patch grid'); end
    fprintf(1, '\n');
    
    %--------------------------------
    %get the connected component (all white regions)
    cctagmap=bwconncomp(tagmap);
    
    % example of contained of cctagmap
    % Connectivity: 8
    % ImageSize: [10234 12174]
    % NumObjects: 23
    % PixelIdxList: {1x23 cell} %1-by-NumObjects cell array where the k-th element in the cell array is a vector containing the linear indices of the pixels in the k-th object.
    
    pp=cctagmap.PixelIdxList;
    h=size(pp);
    ppos=[];
    progress=0;
    for j=1:h(2)
        
        pIdx=pp{1,j};
        [xPosAux yPosAux]=ind2sub(size(tagmap), pIdx(:)); % get the coordenates
        cposAux=0;
        ovlaux=[];
        for k=1:length(xPosAux)
            for m=1:neuronscnt, % go through the annotated neurons   
                ovlp = overlap(xPosAux(k),yPosAux(k),D,D, neurons(m,2),neurons(m,1),neurons(m,4),neurons(m,3));
                if ovlp>=OVERLAP_PIX,
                    ovlaux=[ovlaux;ovlp/(D*D)];
                end
            end
            cposAux=cposAux+1;
            wPosAux(1,cposAux) = max(ovlaux);
        end
        
        %failed idea: try to get the weights from wpos (I tried with
        %differents combinations but the idea is to get the
        %intersection between xpos (all positive pixels) and xPosAux
        %(positive pixels belong to connected component)
        %         [XPOS,idxpos,idxPosAux] = intersect(xpos,xPosAux,'stable');
        %         [YPOS,idypos,idyPosAux] = intersect(ypos,yPosAux,'stable');
        %         [C,itot,icc]=intersect(totalPOS,ccPOS, 'rows','stable');
        %         wPosAux=wpos(itot)';
        %         idx = datasample(totalPOS(itot,:), NPOS, 'Replace', true, 'Weights', wPosAux);
       
        [~,idx] = datasample(xPosAux, NPOS, 'Replace', true,'Weights', wPosAux);
        % go through the selected indexes and extract features from the corresponding patches (x,y,D,D)
     
        for k=1:length(idx),
            % extract positive, save the patch and store to the pdsc/ptag list
            
            patch=imread(imgname,...
                'PixelRegion',{[xPosAux(idx(k)),xPosAux(idx(k))+D-1],[yPosAux(idx(k)),yPosAux(idx(k))+D-1]});
            
            if 0, %IMG_EXPORT,
                pname=sprintf('%s.%d.%d.%d.tif',filearray(i).name,xPosAux(idx(k)),yPosAux(idx(k)),D);
                imwrite(patch, pname);
            end
            
            coord=[coord;str2double(nimg(2)) xPosAux(idx(k)) yPosAux(idx(k))];
            %to save the coordenates in an external file .csv
            auxArray=[str2double(nimg(2)) xPosAux(idx(k)) yPosAux(idx(k))];
            dlmwrite('coordPatches.csv',auxArray,'-append','delimiter',',');
            
            ppos=[ppos; xPosAux(idx(k)) yPosAux(idx(k))];
        
            % extract sift
            patch=single(im2double(patch)*255); % scale to [0,255] and convert to single
            % max val is not 255 but the max encountered within the patch
            [~,d] = vl_sift(patch); % f(not used): 4 x #frames, d: 128 x #frames
            % add it to the patch descriptor list
            pcnt=pcnt+1;
            pdsc{1,pcnt}=d; % concatenate #frames x 128 matrix
            ptag(1,pcnt)=1; % class tag
            
        end
        clear wPosAux xPosAux yPosAux;
        progress1=floor(((j/(h(2)))/0.1));
        if progress1~=progress,
            progress=progress1;
            fprintf(1,'%d%% ', progress*10);
        end
    end
    
    %plot the positive pixels which were choosen
    figure, imshow(tagmap,[]);
    hold on;
    for l=1:length(ppos),
       plot(ppos(l,2),ppos(l,1),'ro');
    end
    export_fig(sprintf('%s.tagmap.%s.%s.pdf',filearray(i).name,'plotPOS',SIGNATURE), '-pdf', '-a2');
    close all;
   
    %save the variables per each mosaic in a .mat-file
    auxname=['m0' num2str(datamat(i))];
    save([auxname '.mat'], 'coord', 'pdsc', 'ptag','D');
    
    fprintf(1,' %6.2f sec\n',toc);
    
    if IMG_EXPORT,
        figure, imshow(imcomplement(img),[]);
        hold on;
        for j=1:size(neurons,1),
            rectangle('Position', neurons(j,:), 'LineWidth', 1, 'EdgeColor', 'y');
        end
        hold on;
        plot(ypos(idx),xpos(idx), 'r.');
        hold on;
        plot(yneg,xneg,'b.');
        
        export_fig(sprintf('%s.trainsmpl.%s.pdf',filearray(i).name,SIGNATURE), '-pdf', '-a2');
        close all;
        
        figure, imshow(tagmap, []);
        export_fig(sprintf('%s.tagmap.%s.pdf',filearray(i).name, SIGNATURE), '-pdf', '-a2');
        close all;
    end
    
    clear img tagmap pdsc ptag coord ; % release memory and clear the variables to be able to reuse them
    pcnt=0;
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