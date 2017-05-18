function [masks,mask_distance,angle_offset,mi0] = dist_seg(images,image_stack,...
    pixel_length,bin_width,scan_mag)

%%
% % mean_im0 = mean2(image_stack);
% % std_im0 = std2(image_stack);
% % thresh = .3*max(image_stack(:));
% % thr_mask = image_stack>thresh;
% % 
% % %apply erosion and then dilation
% % mi0_pre = imfill(thr_mask,'holes');
% % se_erode = strel('disk',1,0);
% % er_mask = imerode(mi0_pre,se_erode);
% % se_dialate = strel('disk',1,0);
% % mi0 = imdilate(er_mask,se_dialate);
% % mi0 = imfill(mi0,'holes');

[mi0] = make_masks(image_stack,1,1.5,10:13,images);
s = regionprops(mi0,'centroid','MajorAxisLength','MinorAxisLength');

if length(s)~=1
    [mi0] = make_masks(image_stack,1,1.25,10:13,images);
    s = regionprops(mi0,'centroid','MajorAxisLength','MinorAxisLength');
end

if length(s)~=1
    [mi0] = make_masks(imgaussfilt(image_stack,2),floor(scan_mag),1.5,10:13);
    if scan_mag==2
        mi0(90:128,:)=0;
        mi0(1:30,:)=0;
        mi0(:,90:128)=0;
        mi0(:,1:30)=0;
    end
    s = regionprops(mi0,'centroid','MajorAxisLength','MinorAxisLength');
end

if length(s)~=1 || ( sum(mi0(:))/16384 > 0.5 && scan_mag~=12.8)
    [mi0] = make_masks(imgaussfilt(image_stack,2),floor(scan_mag),2,10:13);
    if scan_mag==2
        mi0(90:128,:)=0;
        mi0(1:30,:)=0;
        mi0(:,90:128)=0;
        mi0(:,1:30)=0;
    end
    s = regionprops(mi0,'centroid','MajorAxisLength','MinorAxisLength');
end

%[mi0] = make_masks(image_stack,1,scan_mag,0);
%s = regionprops(mi0,'centroid','MajorAxisLength','MinorAxisLength');

if length(s)~=1
    fprintf('multiple regions!!! change threshold');
    dbstop in dist_seg at 42;
    pause;
end

%%
dist_mat = bwdist(bwperim(mi0));
dist_mat = dist_mat.*(~mi0-mi0);
distance_step = bin_width/pixel_length;
mask_distance_end =  min(dist_mat(:))+distance_step:distance_step:max(dist_mat(:));
mask_distance_start = min(dist_mat(:)):distance_step:max(dist_mat(:))-distance_step;
mask_distance = pixel_length * ( mask_distance_start + mask_distance_end ) / 2;

masks = zeros(128,128,length(mask_distance_end));
for k = 1:length(mask_distance_end)
    for i = 1:128
        for j = 1:128
            if mask_distance_start(k) <= dist_mat(i,j) && dist_mat(i,j) < mask_distance_end(k)
                masks(i,j,k) = 1;
            elseif dist_mat(i,j) == mask_distance_end(k) && k == length(mask_distance_end)
                %This condition is for the very furthest pixel
                masks(i,j,k) = 1;      
            end
        end
    end
end

%%
figure(3); clf;
B= imoverlay(mat2gray(image_stack),max(bwperim(mi0),max( masks(:,:,1),...
    masks(:,:,ceil(end/2) ))));
imshow(B,'InitialMagnification','fit');
title('image_stack with 3rd and middle mask');
drawnow;
beep
angle_offset=0;
end











% dist_mat_outward = bwdist(mi0);
% distance_step = bin_width/pixel_length;
% mask_distance_end =  1.5*distance_step:distance_step:max(dist_mat_outward(:));
% mask_distance_start = distance_step/2:distance_step:max(dist_mat_outward(:))-distance_step;
% mask_distance_outward = pixel_length * ( mask_distance_start + mask_distance_end ) / 2;
% %%
% masks_outward = zeros(128,128,length(mask_distance_end));
% for k = 1:length(mask_distance_end)
%     for i = 1:128
%         for j = 1:128
%             if mask_distance_start(k) <= dist_mat_outward(i,j) && dist_mat_outward(i,j) < mask_distance_end(k)
%                 masks_outward(i,j,k) = 1;
%             elseif dist_mat_outward(i,j) == mask_distance_end(k) && k == length(mask_distance_end)
%                 %This condition is for the very furthest pixel
%                 masks_outward(i,j,k) = 1;      
%             end
%         end
%     end
% end
% %removes first masks (filled in region, mi0) from masks outward
% %masks_outward = masks_outward(:,:,2:end);
% 
% dist_mat_inward = bwdist(~max(mi0,bwperim(mi0)));
% mask_distance_end =  1.5*distance_step:distance_step:max(dist_mat_inward(:));
% mask_distance_start = distance_step/2:distance_step:max(dist_mat_inward(:))-distance_step;
% mask_distance_inward = fliplr(-pixel_length *( mask_distance_start + mask_distance_end ) / 2);
% masks_inward = zeros(128,128,length(mask_distance_end));
% for k = 1:length(mask_distance_end)
%     k_ind = length(mask_distance_end)-k+1;
%     for i = 1:128
%         for j = 1:128
%             if mask_distance_start(k) <= dist_mat_inward(i,j) && dist_mat_inward(i,j) < mask_distance_end(k)
%                 masks_inward(i,j,k_ind) = 1;
% %             elseif dist_mat_inward(i,j) == mask_distance_end(k) && k == length(mask_distance_end)
% %                 %This condition is for the very furthest pixel
% %                 masks_inward(i,j,k_ind) = 1;      
%             end
%         end
%     end
% end
% %removes first masks (filled in region, ~mi0) from masks outward
% %masks_inward = masks_inward(:,:,2:end);
% 
% dist_mat_center = bwdist(bwperim(mi0));
% mask_distance_end =  distance_step/2;
% mask_distance_start = 0;
% masks_center = zeros(128,128,1);
% for i = 1:128
%     for j = 1:128
%         if mask_distance_start <= dist_mat_center(i,j) && dist_mat_center(i,j) < mask_distance_end
%             masks_center(i,j,1) = 1;
%         end
%     end
% end



%masks = cat(3,masks_inward,masks_center,masks_outward);
%mask_distance = [mask_distance_inward,0, mask_distance_outward];
