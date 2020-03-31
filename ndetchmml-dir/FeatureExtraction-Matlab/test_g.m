function test_g(S, NLOCS, OPOS, TRAINDATA_PATH)
% clc;
% clear all;
% close all;
%--------------------------------------
% take annotated mosaics from the TEST_DIR (m01.tif, m01.log;...)
% feature extraction - extract square test patches in each mosaic:
% coordinate (patch) sampling - it is important to capture the mosaic- sampling every pixel
% is computationally heavy
% (option 1): sample a regular grid (step=S*D along x and y) that captures the whole mosaic
% (option 2): uniform random sample coordinate along x and y
% (option 3): random sample coordinate with replacement using intensity as importance
% sampling function
% sample the patches, extract sift descriptors and codeword histograms,
% apply the trained classifiers - visualize and evaluate
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
%load('./experiments/traindata.mat'); % path to traindata.mat
load(TRAINDATA_PATH);%('./TRAIN,D=500,S=0.5,OPOS=0.5,NW=50/traindata.mat');
[~, TRAINDATA_NAME, ~]=fileparts(TRAINDATA_PATH);
%--------------------------------------
% D, S, OPOS, NW, centroids loaded from traindata.mat
NW=size(centroids,2);

% define PATCH_SAMPLING_METHOD 0:grid, 1:random uniform 2(experimental):random image weighted
if S>0,
    PATCH_SAMPLING_METHOD=0; % grid, S*D step along x,y, NLOCS is irrelevant
    SIGNATURE=sprintf('S.OPOS_%3.1f_%3.1f_USING[%s]', S, OPOS,TRAINDATA_NAME);
else
    PATCH_SAMPLING_METHOD=1; % random uniform, NLOCS used
    SIGNATURE=sprintf('NLOCS.OPOS_%d_%3.1f_USING[%s]', NLOCS, OPOS,TRAINDATA_NAME);
end

disp(SIGNATURE);

% S, NLOCS, OPOS are read as arguments

OVERLAP_PIX=OPOS*D*D; % criteria to consider location as positive
IMG_EXT='tif'; % mosaic extension
IMG_EXPORT=1; % save midresults

% check if directory and files exist
TEST_DIR='C:\Users\gamata\Pictures\mosaics\test'; % name of the subdir with test tifs/logs m01.tif,m01.log,...
if isdir(TEST_DIR) == 0, error('TEST directory does not exist'); end
%--------------------------------------
filearray = dir([TEST_DIR filesep '*.' IMG_EXT]); % get all files in the directory
NumImages = size(filearray,1); % get the number of images
fprintf(1,'%d images found in %s\n', NumImages, TEST_DIR);
if NumImages < 0, error('No image in the directory'); end

% initially allocate 2000 patches per train image
phist=zeros(2000*NumImages, NW); % patch histograms
ptag=zeros(1,2000*NumImages);    % patch class tag (0=BACKGROUND,1=NEURON)
pcnt=0;                          % patch count

for i=1:NumImages,
    
    imgname = [TEST_DIR filesep filearray(i).name]; % read test mosaic
    img=imread(imgname);
    
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
    tagmap=zeros(size(img)); % to discard positives while going through the grid
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
                    end
                end
            end
        end
    end
    
    fprintf(1,'%s, [%d, %d] vl_sift():\n', imgname, size(img,1), size(img,2));
    
    % Xtest, Ytest
    if PATCH_SAMPLING_METHOD==0,
        [X, Y]=ndgrid(randi([1 round(0.3*D)]):round(S*D):size(img,1)-D+1, randi([1 round(0.3*D)]):round(S*D):size(img,2)-D+1);
        Xtest=X(:);
        Ytest=Y(:);
        clear X Y;
    else if PATCH_SAMPLING_METHOD==1,
            % random uniform
            Xtest=randi([randi([1 round(0.1*D)]) size(img,1)-D+1], 1, NLOCS);
            Ytest=randi([randi([1 round(0.1*D)]) size(img,2)-D+1], 1, NLOCS);
        else if PATCH_SAMPLING_METHOD==2, % (experimental)
                % random image weighted
                [X,Y]   = ndgrid(1:2:size(img,1)-D+1, 1:2:size(img,2)-D+1);
                w=img(1:2:size(img,1)-D+1, 1:2:size(img,2)-D+1).^2;
                idx = datasample(1:length(X(:)), NLOCS, 'Replace', true, 'Weights', w(:));
                
                Xtest=X(idx);
                Ytest=Y(idx);
                clear X Y w;
            end
        end
    end
    
    ppi_all=length(Xtest); % nr patches per image
    ppi_cnt=0; % patch per image count (to plot progress)

    progress=0; 
    
    beg_idx=pcnt;
    
    tic;
    for j=1:length(Xtest),    
        patch=imread(imgname,'PixelRegion',{[Xtest(j),Xtest(j)+D-1],[Ytest(j),Ytest(j)+D-1]});
            
        ppi_cnt = ppi_cnt+1;
        progress1=floor(((ppi_cnt/ppi_all)/0.1));
        if progress1~=progress,
        	progress=progress1;
            fprintf(1,'%d%% ', progress*10);
        end
            
        % extract sift
        patch=single(im2double(patch)*255); % scale to [0,255] and convert to single
        % max val is not 255 but the max encountered within the patch
        [~,d] = vl_sift(patch); % f(not used): 4 x #frames, d: 128 x #frames
        [~, k] = min(vl_alldist(double(d), centroids), [], 2); % assign centroid to each patch descriptor
            
        pcnt=pcnt+1;
        phist(pcnt,:)=hist(k,1:NW);
        ptag(1,pcnt)=tagmap(Xtest(j),Ytest(j)); % class tag
            
    end
    fprintf(1,' %6.2f sec\n',toc);
    end_idx=pcnt;
    
    % done extracting features for current mosaic - plot detections
    A=dataset(phist(beg_idx+1:end_idx,:),ptag(1,beg_idx+1:end_idx)');% take those from the last image
    Ac=A*W*labeld; % get the classification results for this train image (vizualization)
    fprintf(1, '\n');
    %---------------manually---------------------------
    %recall, precision,... for each classifier
    
    tt=ptag(1,beg_idx+1:end_idx)';%real tags
   
    TPTN=[];
    TP=[];
    TN=[];
    FP=[];
    FN=[];
    PredP=[];
    recall=[];
    precision=[];
    accuracy=[];
    
    %total and P are the same numbers for each classifier 
   
    P=sum(tt); %real positive
    
    for k=1:length(Ac),
       clas=Ac{k};
       
       %confusion matrix (some can be calculated based on others, but to be
       %sure that everythings are working right)
       TP=[TP;sum((tt==clas)&(tt==1))];%true positive
       TN=[TN;sum((tt==clas)&(tt==0))];%true negative
       FP=[FP;sum((tt~=clas)&(tt==0))];%false positive
       FN=[FN;sum((tt~=clas)&(tt==1))];%false negative
       PredP=[PredP;sum(clas)];%total positive predicted
       TPTN=[TPTN;sum(tt==clas)];%TP plus TN
       %evaluation
       recall=[recall;TP(length(TP))/P(length(P))];%TP/P
       precision=[precision;TP(length(TP))/PredP(length(PredP))];%TP/(TP+FP) = TP/PredP
       accuracy=[accuracy;(TP(length(TP))+TN(length(TN)))/(length(tt))];%(TP+TN)/nºpatches
      
        
    end
    AcAux=Ac';
    cnamesAux=cnames';
    save('confusion_matrix.mat','tt','AcAux' ,'cnamesAux', 'TPTN', 'TP', 'TN', 'P','FP','FN','PredP','recall','precision','accuracy');
    %--------------------------------------------------
    if IMG_EXPORT,
        for k=1:length(Ac), % for each classifier 
            % output an image with test patch sampling and detections overlayed
            figure, imshow(imcomplement(img),[]);
            hold on;
            plot(Ytest,Xtest, 'y.');
            hold on;
            % go through the classifier output, add detection rectangles
            for l=1:length(Xtest),
                if Ac{k}(l)==1,
                    rectangle('Position', [Ytest(l) Xtest(l) D D], 'LineWidth', 1, 'EdgeColor', 'r');
                end
            end
            export_fig(sprintf('%s.%s.%s.pdf',filearray(i).name,cnames{k},SIGNATURE), '-pdf', '-a2');
            close all;
            
            figure, imshow(tagmap,[]);
            hold on;
            plot(Ytest,Xtest, 'y.');
            hold on;
            % go through the classifier output, add detection rectangles
            for l=1:length(Xtest),
                if Ac{k}(l)==1,
                    %rectangle('Position', [Ytest(l) Xtest(l) D D], 'LineWidth', 1, 'EdgeColor', 'r');
                    plot(Ytest(l),Xtest(l),'ro');
                end
            end
            export_fig(sprintf('%s.tagmap.%s.%s.pdf',filearray(i).name,cnames{k},SIGNATURE), '-pdf', '-a2');
            close all;
        end
    end
    
    clear img tagmap;
        
end

fprintf(1,'done extracting features, [%d test objects]\n', pcnt);

% crop unused allocated patch descriptors and tags if existent
phist=phist(1:pcnt,:);
ptag=ptag(1,1:pcnt);

% evaluation
A=dataset(phist, ptag');
%AWl=A*W*labeld; % get the classification results for this train image (vizualization)
 
% precision
Pre=A*W*testc([],'precision');
figure, bar(cell2mat(Pre)); set(gca,'xticklabel',cnames)
export_fig(sprintf('p.%s.pdf',SIGNATURE), '-pdf', '-a2', '-native', '-transparent');
close all;

% recall
Rec=A*W*testc([],'sensitivity');
figure, bar(cell2mat(Rec)); set(gca,'xticklabel',cnames)
export_fig(sprintf('r.%s.pdf',SIGNATURE), '-pdf', '-a2', '-native', '-transparent');
close all;

% classification error based on error counting,  weighted by the class priors
Err1=A*W*testc([],'crisp');
figure, bar(cell2mat(Err1)); set(gca,'xticklabel',cnames)
export_fig(sprintf('err1.%s.pdf',SIGNATURE), '-pdf', '-a2', '-native', '-transparent');
close all;

% classification error based on soft error  summation, i.e. a sum of the absolute difference between  classifier output and target, weighted by class priors
Err2=A*W*testc([],'soft');
figure, bar(cell2mat(Err2)); set(gca,'xticklabel',cnames)
export_fig(sprintf('err2.%s.pdf',SIGNATURE), '-pdf', '-a2', '-native', '-transparent');
close all;

% Area under the ROC curve (this is an error and not a  performance!)
Auc=A*W*testc([],'auc');
figure, bar(cell2mat(Auc)); set(gca,'xticklabel',cnames)
export_fig(sprintf('auc.%s.pdf',SIGNATURE), '-pdf', '-a2', '-native', '-transparent');
close all;

testfeats=[phist ptag'];
save(sprintf('testdata.%s.mat',SIGNATURE), 'testfeats', 'Pre', 'Rec', 'Err1', 'Err2', 'Auc', 'D', 'S', 'NLOCS', 'OPOS');    
csvwrite(sprintf('testfeats.%s.csv', SIGNATURE), testfeats);

% if PATCH_SAMPLING_METHOD==0,
%     save(sprintf('testdata.%s.mat',SIGNATURE), 'testfeats', 'Pre', 'Rec', 'Err1', 'Err2', 'Auc', 'D', 'S', 'NLOCS', 'OPOS');
%     csvwrite(sprintf('testfeats.%s.csv', SIGNATURE), testfeats);
% else
%     save(sprintf('testdata.%s.mat',SIGNATURE), 'testfeats', 'Pre', 'Rec', 'Err1', 'Err2', 'Auc', 'D', 'NLOCS', 'OPOS');
%     csvwrite(sprintf('testfeats.D.NLOCS.OPOS.NW_%d_%d_%3.1f_%d.csv', D, NLOCS, OPOS, NW), testfeats);
% end

fprintf(1,'finished.\n');

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

