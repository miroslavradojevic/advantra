MOSAIC_DIR = '/Users/miroslav/nmachinel.test/data'; % dir with mXY.tif
D = 500;
S = 1.0;
% run descriptor_mosaic for each mXY.tif
SIGNATURE=sprintf('D_S_%d_%3.1f', D, S);
filearray = dir([MOSAIC_DIR filesep '*.tif']); % get all files in the directory
NumImages = size(filearray,1); % get the number of images

for i=1:NumImages,
    imgname = [MOSAIC_DIR filesep filearray(i).name]; % read train mosaic
    [fdir,fname,fext] = fileparts(imgname);
%     fprintf(1,'%s\n',fdir);
%     fprintf(1,'%s_dmap_%s.mat\n', fname, SIGNATURE)
    fprintf(1,'%s\n',imgname);
    dmap = descriptor_mosaic(imgname, D, S);
    save([fdir filesep sprintf('%s_dmap_%s.mat', fname, SIGNATURE)],'dmap');
end


% indeed this can be executed for each mosaics and called in parallel
% trhough terminal matlavb call or even opening several matlab gui
% instances
