FOV_pth = 'na'; FOV_fn = 'none'; tsm = 0; reach = 0; tmini=1; tmaxi = 4096;
scan_mag = 2;
data_path = '/Users/bryankaye/Documents/MATLAB/data/2017-02-23/';
data_name ='2X_R1_S1';
[tsm,dataname_cell] = find_filenames(data_path,data_name,tsm,1);
num_images = length(dataname_cell);

for k = 1:num_images
    [int_image0{k},FLIMages{k}] = spc2FLIMdata(tmini,tmaxi,dataname_cell{k},data_path,FOV_fn,FOV_pth,tsm,reach);
end
rot_method = 'bicubic';
image_stack = int_image0{k};
pixel_length = (440 / scan_mag) / 128;
[image_stack,registration_vectors] = register_images(int_image0,...
    scan_mag,rot_method);
%%
n_masks = 40;
mean_im0 = mean2(image_stack);
std_im0 = std2(image_stack);
thresh = .2*max(image_stack(:));
thr_mask = image_stack>thresh;

%apply erosion and then dilation
mi0_pre = imfill(thr_mask,'holes');
se_erode = strel('disk',1,0);
er_mask = imerode(mi0_pre,se_erode);
se_dialate = strel('disk',1,0);
mi0 = imdilate(er_mask,se_dialate);
mi0_props = regionprops(mi0,'centroid','MajorAxisLength',...
    'MinorAxisLength','Orientation');
centroid(1) = mi0_props.Centroid(1)-65;
centroid(2) = mi0_props.Centroid(2)-65;

mi = mi0;
figure(3); clf; imshow(mat2gray(thr_mask.*image_stack),'InitialMagnification','fit');
drawnow;

dist_mat = bwdist(thr_mask);

distance_step = max(dist_mat(:))/n_masks;
mask_distances_end =  distance_step:distance_step:max(dist_mat(:));
mask_distances_start = 0:distance_step:max(dist_mat(:))-distance_step;
mask_distances = ( mask_distances_start + mask_distances_end ) / 2;
%%
tic
for k = 1:length(mask_distances_end)
    masks{k} = zeros(128,128);
    for i = 1:128
        for j = 1:128
            if mask_distances_start(k) <= dist_mat(i,j) && dist_mat(i,j) < mask_distances_end(k)
                masks{k}(i,j) = 1;
                %add else for last data point elseif distmat(i,j)
                %==mask_distance_end(end) && k = kmax. Currently we lose very furthest pixel       
            end
        end
    end
end
toc
%%
h3 = figure(3); clf;
B= imoverlay(mat2gray(image_stack),max(masks{1},masks{20}));
imshow(B,'InitialMagnification','fit');
%imagesc(int_final)
%axis equal
%hold on;
title('aligned registered image and first and last mask');
angle_offset=0;