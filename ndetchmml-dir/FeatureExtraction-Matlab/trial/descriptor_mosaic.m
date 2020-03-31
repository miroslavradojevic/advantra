function dmap = descriptor_mosaic(mosaic_path, D, S)
% input uint8 image mosaic
% D patch scale
% S overlap (def. sampling density) 
% output set of desriptors assigned to each (x,y) location of the mosaic
VLFEAT_DIR='/Users/miroslav/vlfeat-0.9.20'; % vlfeat library path
run(fullfile(VLFEAT_DIR,'toolbox/vl_setup')); % link vlfeat library
% vl_version verbose % check if library is ok

m=imread(mosaic_path);
dmap=cell(size(m)); % convenient to store in 2d cell for faster readout later in descriptor_patch() 

Xgrid = 1 : round(S*D) : size(m,1)-D+1;
Ygrid = 1 : round(S*D) : size(m,2)-D+1;
ppi_all=length(Xgrid)*length(Ygrid); % nr patches per image

clear m % release memory

% go through the mosaic
ppi_cnt=0; % patch per image count (to plot progress)
    
progress=0; 
tic; % measure time
for x = Xgrid,
    for y = Ygrid,
        patch=imread(mosaic_path,'PixelRegion',{[x,x+D-1],[y,y+D-1]});
            
%         imshow(imcomplement(patch));
            
            ppi_cnt = ppi_cnt+1;
            progress1=floor(((ppi_cnt/ppi_all)/0.05));
            if progress1~=progress,
                progress=progress1;
                fprintf(1,'%d%% ', progress*5);
            end
            
            patch=single(im2double(patch)*255); % scale to [0,255] and convert to single

            [frm,dsc] = vl_sift(patch); % f(not used): 4 x #frames, d: 128 x #frames
            
            for i = 1 : size(frm,2),
                xF = floor(frm(1,i));
                yF = floor(frm(2,i));
                dmap{xF, yF}(size(dmap{xF, yF},1)+1, : ) = dsc(:,i)';
            end  
    end % Ygrid
end % Xgrid
fprintf(1,' %6.2f sec\n',toc);
