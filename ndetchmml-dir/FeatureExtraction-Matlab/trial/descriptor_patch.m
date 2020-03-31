function d = descriptor_patch(dmap, x, y, W, H)
% does what vl_sift was - retrieves descriptors
% within some rectangle
% this one plugged into train.m at the place where vl_sift command is
d= zeros(1000,128); % parch descriptor
cnt = 0;
for x1 = x : x+W,
    for y1 = y : y+H,
        if size(dmap{x1, y1},1)>0,
            d(cnt+1:cnt+size(dmap{x1,y1},1),:) = dmap{x1,y1};
            cnt = cnt + size(dmap{x1,y1},1);
        end       
    end
end
d=d(1:cnt,:); % crop
fprintf(1,'%d descriptors found\n',cnt);