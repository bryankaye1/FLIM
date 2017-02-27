function [masks,mask_distance,angle_offset] = dist_seg(image_stack,pixel_length)

%%
n_masks = 40;
mean_im0 = mean2(image_stack);
std_im0 = std2(image_stack);
thresh = .3*max(image_stack(:));
thr_mask = image_stack>thresh;

%apply erosion and then dilation
mi0_pre = imfill(thr_mask,'holes');
se_erode = strel('disk',1,0);
er_mask = imerode(mi0_pre,se_erode);
se_dialate = strel('disk',1,0);
mi0 = imdilate(er_mask,se_dialate);
figure(5); clf; imshow(mi0);
figure(6); clf; imshow(thr_mask);

mi0_props = regionprops(mi0,'centroid','MajorAxisLength',...
    'MinorAxisLength','Orientation');
if length(mi0_props)>1
    input('multiple regions!!! change threshold');
end

figure(5); clf; imshow(mat2gray(mi0.*image_stack),'InitialMagnification','fit');

dist_mat = bwdist(mi0);
distance_step = max(dist_mat(:))/n_masks;
mask_distance_end =  distance_step:distance_step:max(dist_mat(:));
mask_distance_start = 0:distance_step:max(dist_mat(:))-distance_step;
mask_distance = pixel_length*( mask_distance_start + mask_distance_end ) / 2;
%%
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
B= imoverlay(mat2gray(image_stack),max( masks(:,:,3),masks(:,:,ceil(end/2) )));
imshow(B,'InitialMagnification','fit');
%imagesc(int_final)
%axis equal
%hold on;
title('image_stack with 3rd and middle mask');
drawnow;
beep
angle_offset=0;
end
