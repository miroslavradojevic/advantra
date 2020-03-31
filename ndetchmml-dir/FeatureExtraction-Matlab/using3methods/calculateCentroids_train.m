function calculateCentroids_train(NW)
 
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
%---------------------------
IMG_EXPORT=1; %save midresults
%.D.S.OPOS.NW.NPOS_%d_%3.1f_%3.1f_%d_%d
D=500; % rectangular patch size: width/height
S=0.5; % defines patch sampling grid step relative to patch size D step=S*D
OPOS=0.5; % overlap to be considered as positive
NW=60; % number of clusters, codewords to form the dictionary with
NPOS=33; % number of patches for each patch marked like neuron by Miguel.
SIGNATURE=sprintf('D.S.OPOS.NW.NPOS_%d_%3.1f_%3.1f_%d_%d', D, S, OPOS, NW, NPOS); %.D.S.OPOS.NW.NPOS_500_0.5_0.5_75_33;
%-----------------------------
%number of mosaic to train
datamat=[1 2 3 5 6 7];
NumImages=length(datamat);

tstart=tic;

%put all descriptors (from all patches, from all mosaics) in the same cell
auxMat=matfile(['m0' num2str(datamat(1)) '.mat']);
pdsc=auxMat.pdsc;
for i=2:NumImages
    auxMat=matfile(['m0' num2str(datamat(i)) '.mat']);
    pdsc=[pdsc auxMat.pdsc];
end

% calculate centroids using kmeans find the best dictionary or codebook to vector quantize the data.
fprintf(1, 'vl_kmeans()...');
tic;
[centroids, ~] = vl_kmeans(double(cell2mat(pdsc)), NW);
fprintf(1, '\n\t%6.2f sec\n',toc);

% use the centroids to calculate the histogram for each patch
% go through the patches again and express each wrt dictionary
trainfeats=zeros(length(pdsc),NW);
for i=1:length(pdsc),
    [~, k] = min(vl_alldist(double(pdsc{i}), centroids), [], 2); % assign centroid to each patch descriptor
    trainfeats(i,:)=hist(k,1:NW); % distribution of the centroids in current image patch
end

%put all tags (from all patches, from all mosaics) in the same cell
auxMat=matfile(['m0' num2str(datamat(1)) '.mat']);
ptag=auxMat.ptag;
for i=2:NumImages
    auxMat=matfile(['m0' num2str(datamat(i)) '.mat']);
    ptag=[ptag auxMat.ptag];
end

trainfeats=[trainfeats ptag']; % intentionally left out (xy from the feature set here as those are important to be recorded in the test stage)
% also training a classifier at this stage using x and y would make a mess
% as they are not corerct features on the other hand it is useful (and necessary to have xy added when calculating features in the test stage to interpret the results, visualize, regress)

if IMG_EXPORT,
    % take up to 5 random positives/negatives and visualize their codeword histograms
    positives=find(trainfeats(:,end)==1);
    perm = randperm(length(positives));
    sel = perm(1:min(5,length(positives)));
    
    figure,
    for i=1:length(sel),
        subplot(5,1,i);
        bar(trainfeats(positives(sel(i)),1:end-1));
    end
    export_fig(sprintf('pos5.%s.pdf',SIGNATURE), '-pdf', '-a2', '-native', '-transparent');
    close all;
    
    negatives=find(trainfeats(:,end)==0);
    perm = randperm(length(negatives));
    sel = perm(1:min(5,length(negatives)));
    
    figure,
    for i=1:length(sel),
        subplot(5,1,i);
        bar(trainfeats(negatives(sel(i)),1:end-1));
    end
    export_fig(sprintf('neg5.%s.pdf',SIGNATURE), '-pdf', '-a2', '-native', '-transparent');
    close all;
    
end

fprintf(1,'done extracting features, %d[%d positives]\n', length(ptag), size(find(ptag==1),2));

% save ouput, append features and centroids to train.mat
% features: matrix with columns: h[1] h[2]... h[NW] classtag
% centroids: 128 x NW matrix (with sift descriptors along columns)
aux=sprintf('traindata.%s.mat',SIGNATURE);
save(aux,'trainfeats','centroids','D', 'S', 'OPOS');%, 'W', 'cnames', 'e10');
csvwrite(sprintf('trainfeats.%s.csv', SIGNATURE), trainfeats);
csvwrite(sprintf('centroids.%s.csv', SIGNATURE), centroids);
fprintf(1, 'finished.\n');
tend2=toc(tstart);
fprintf(1, '\n\t%6.2f sec\n',tend2);
