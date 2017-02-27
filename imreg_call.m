%This file is used to test image segmentation
%Written by Olivia Stiehl and Bryan Kaye
clear;
sim_im = 0;
%num_images = 20;%length(dataname_cell);
rot_method = 'bicubic';
if ~sim_im
    FOV_pth = 'na';
    FOV_fn = 'none';
    tsm = 0; reach = 0;
    tmini=1; tmaxi = 4096;
    data_path = '/Users/bryankaye/Documents/MATLAB/data/2017-01-26/';
    data_name ='Donor_r1';
    [tsm,dataname_cell] = find_filenames(data_path,data_name,tsm,1);
    num_images = length(dataname_cell);
    int_image0 = cell(1,num_images);
    for k = 1:num_images
        [int_image0{k}] = spc2image(tmini,tmaxi,dataname_cell{k},data_path,FOV_fn,FOV_pth,tsm,reach);
        [FLIMage_temp] = sdt_to_vector(data_path,dataname_cell{k});
        FLIMages{k} = shiftdim(FLIMage_temp,1);
        %need to reshape FLIMage ld(tmini:tmaxi,:,:);
    end    
else
    data_path = '/Users/bryankaye/Documents/MATLAB/data/sim_images/';
    base_name = 'sim_image';
%    int_image0{1} = imread([data_path,'sim_image1.tif']);
%    int_image0{2} = imread([data_path,'sim_image2.tif']);
     namelist_str = dir([data_path,base_name,'*']);
     for dc_ind = 1:length(namelist_str)
        dataname_cell{dc_ind} = namelist_str(dc_ind).name;
        int_image0{dc_ind} = imread([data_path,dataname_cell{dc_ind}]);
     end
     num_images = length(namelist_str);
end
%%
[image_stack,registration_vectors] = register_images(int_image0,rot_method);
[int_masks,mask_angle_rot] = ellips_seg(image_stack,rot_method);
[FLIM_stack] = transform_FLIMage(FLIMages,registration_vectors,...
    mask_angle_rot,rot_method);
%%
for k = 1:size(int_masks,3)
  FLIM_masks{k} = repmat(int_masks(:,:,k),1,1,tmaxi-tmini+1);
  FLIMseg{k} = FLIM_masks{k}.*  FLIM_stack;
  FLIMgroups{k} = squeeze(sum(sum(FLIMseg{k}))); %3rd dimension is FLIM histogram
  
  intseg{k} = int_masks(:,:,k).* imrotate(image_stack,-mask_angle_rot,rot_method,'crop');
  intgroups{k} = squeeze(sum(sum(intseg{k}))); %3rd dimension is FLIM histogram
end


%%

for k =1:size(int_masks,3)
    plot_int(k) = intgroups{k};
    plot_FLIM_int(k) = sum(FLIMgroups{k});  
end

figure(4); clf; plot(plot_int,plot_FLIM_int, 'bo');
title('Photon counts match between intensity image and FLIMage?');
xlabel('intensity image photons per mask');
ylabel('FLIMage photons per mask');
%%
% figure(20); clf; 
% for i =1:num_images
%     subplot(1,num_images+1,i); 
%     imshow(mat2gray(int_image0{i}));
% end
% subplot(1,num_images+1,i+1);
% imshow(mat2gray(image_stack));



