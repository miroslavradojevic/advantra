function calculateHistograms_test(TRAINDATA_PATH) %S, OPOS,NLOCS
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
datamat=[1 2 3 4 5 6 7 8];
NumImages=length(datamat);
D=500; % rectangular patch size: width/height
S=0.5; % defines patch sampling grid step relative to patch size D step=S*D
OPOS=0.5; % overlap to be considered as positive
NW=60; % number of clusters, codewords to form the dictionary with
NPOS=33; % number of patches for each patch marked like neuron by Miguel.
SIGNATURE=sprintf('D.S.OPOS.NW.NPOS_%d_%3.1f_%3.1f_%d_%d', D, S, OPOS, NW, NPOS); %.D.S.OPOS.NW.NPOS_500_0.5_0.5_75_33;

%--------------------------------------
%load('./experiments/traindata.mat'); % path to traindata.mat
load(TRAINDATA_PATH);%('./TRAIN,D=500,S=0.5,OPOS=0.5,NW=50/traindata.mat');
[~, TRAINDATA_NAME, ~]=fileparts(TRAINDATA_PATH);
%--------------------------------------
% D, S, OPOS, NW, centroids loaded from traindata.mat
NW=size(centroids,2);

disp(SIGNATURE);

% check if directory and files exist
% TEST_DIR='C:\Users\gamata\Pictures\mosaics\test'; % name of the subdir with test tifs/logs m01.tif,m01.log,...
% if isdir(TEST_DIR) == 0, error('TEST directory does not exist'); end
%--------------------------------------
% filearray = dir([TEST_DIR filesep '*.' IMG_EXT]); % get all files in the directory
% NumImages = size(filearray,1); % get the number of images
% fprintf(1,'%d images found in %s\n', NumImages, TEST_DIR);
if NumImages < 0, error('No image in the directory'); end

% initially allocate 2000 patches per train image
phist=zeros(2000*NumImages, NW); % patch histograms
pcnt=0;                          % patch count

for i=1:NumImages,
    
        tic;
        pcnt=0; 
        auxMat=matfile(['m0' num2str(datamat(i)) '.mat']);
        pdsc=auxMat.pdsc;
        ptag=auxMat.ptag;
        
        for k=1:length(pdsc),
            auxD=pdsc(k);
            d=cell2mat(auxD);
            
            [~, k] = min(vl_alldist(double(d), centroids), [], 2); % assign centroid to each patch descriptor
            
            pcnt=pcnt+1;
             phist(pcnt,:)=hist(k,1:NW);
        end
    fprintf(1,' %6.2f sec\n',toc);
    
    save(['m0' num2str(datamat(i)) '_hist.mat'], 'phist', 'ptag');
    dlmwrite(['m0' num2str(datamat(i)) '_hist.csv'],phist,'-append','delimiter',',');
    clear phist;
        
end

fprintf(1,'done extracting features, [%d test objects]\n', pcnt);
