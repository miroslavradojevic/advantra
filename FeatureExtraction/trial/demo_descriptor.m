function demo_descriptor(mosaic_path, D)
% input uint8 image mosaic
% D patch scale
% output set of desriptors assigned to each (x,y) location of the mosaic
N=10; % number of demos
VLFEAT_DIR='/Users/miroslav/vlfeat-0.9.20'; % vlfeat library path
run(fullfile(VLFEAT_DIR,'toolbox/vl_setup')); % link vlfeat library
% vl_version verbose % check if library is ok

m=imread(mosaic_path);
% dmap=cell(size(m)); % convenient to store in 2d cell for faster readout later in descriptor_patch() 

% Xgrid = 1 : round(S*D) : size(m,1)-D+1;
% Ygrid = 1 : round(S*D) : size(m,2)-D+1;
% ppi_all=length(Xgrid)*length(Ygrid); % nr patches per image

x = randi([1 size(m,1)],1,N);
y = randi([1 size(m,2)],1,N);

clear m % release memory

% go through the mosaic
ppi_cnt=0;
    
% progress=0; 
% tic; % measure time
% for x = Xgrid,
%     for y = Ygrid,



for n=1:length(x),
        patch=imread(mosaic_path,'PixelRegion',{[x(n),x(n)+D-1],[y(n),y(n)+D-1]});
        
        figure;
        imshow(imcomplement(patch));
        
ppi_cnt = ppi_cnt+1;

%         progress1=floor(((ppi_cnt/ppi_all)/0.05));
%         if progress1~=progress,
%             progress=progress1;
%             fprintf(1,'%d%% ', progress*5);
%         end
            
            patch=single(im2double(patch)*255); % scale to [0,255] and convert to single

            tic;
            [frm,dsc] = vl_sift(patch); % f(not used): 4 x #frames, d: 128 x #frames
            toc
%             perm = randperm(size(frm,2)) ;
%             sel = perm(1:50);
        h1 = vl_plotframe(frm); % (:,sel)
        h2 = vl_plotframe(frm); % (:,sel) 
        set(h1,'color','k','linewidth',3);
        set(h2,'color','y','linewidth',2);
            
%             if ppi_cnt>=N,
%                 break
%             end
%             for i = 1 : size(frm,2),
%                 xF = floor(frm(1,i));
%                 yF = floor(frm(2,i));
%                 dmap{xF, yF}(size(dmap{xF, yF},1)+1, : ) = dsc(:,i)';
%             end  
%     end % Ygrid
% end % Xgrid
end
% fprintf(1,' %6.2f sec\n',toc);

