%This file is used to test image segmentation
%Written by Olivia Stiehl and Bryan Kaye

function [FLIMgroups,int_groups,reg_seg_plots] = spindle_area_reg_seg(data_path,data_name,...
    tmini,tmaxi,FOV_fn,FOV_pth,mask_type,scan_mag,pixel_bin_width)


sim_im = 0;
rot_method = 'bilinear'; %just did bilinear
if ~sim_im
    %FOV_pth = 'na'; FOV_fn = 'none'; mask_type = 'ellipsoid'; scan_mag = 8;
    %tmini=1; tmaxi = 4096; data_path = '/Users/bryankaye/Documents/MATLAB/data/2017-01-26/';
    %data_name ='Donor_r1'; 
    [tsm,dataname_cell] = find_filenames(data_path,data_name,0,1);
    num_images = length(dataname_cell);

    for k = 1:num_images     
        [int_image0{k},FLIMages{k}] = spc2FLIMdata(tmini,tmaxi,dataname_cell{k},data_path,FOV_fn,FOV_pth,tsm,reach);
    end

else
    data_path = '/Users/bryankaye/Documents/MATLAB/data/sim_images/';
    base_name = 'sim_image';
    namelist_str = dir([data_path,base_name,'*']);
    for dc_ind = 1:length(namelist_str)
        dataname_cell{dc_ind} = namelist_str(dc_ind).name;
        int_image0{dc_ind} = imread([data_path,dataname_cell{dc_ind}]);
    end
    num_images = length(namelist_str);
end
%%
pixel_length = (440 / scan_mag) / 128;
short_axis_len = 10 ./ ((440./scan_mag) / 128); %minimum spindle short axis length in pixels 
max_conv_dist = floor(short_axis_len*.7);
[image_stack,registration_vectors] = register_images(int_image0,...
    scan_mag,rot_method);

if strcmp(mask_type,'ellipsoid')
[int_masks,mask_distance,mask_angle_rot] = ellips_seg(image_stack,...
    pixel_length,scan_mag,rot_method);
elseif strcmp(mask_type, 'edge_distance')
    [int_masks,mask_distance,mask_angle_rot] = dist_seg(int_image0,image_stack,...
        pixel_length,pixel_bin_width,scan_mag);
else
    input('ERROR: NO MASK TYPE SELECTED'); 
end

[FLIM_stack] = transform_FLIMage(FLIMages,registration_vectors,...
    mask_angle_rot,rot_method);
%%
%clear FLIM_masks FLIMseg FLIMgroups plot_FLIM_int intseg int_seg_temp int_groups
for k = 1:size(int_masks,3)
    FLIM_masks = repmat(int_masks(:,:,k),1,1,tmaxi-tmini+1);
    FLIMseg = FLIM_masks .* FLIM_stack;
    FLIMgroups(k,:) = squeeze(sum(sum(FLIMseg))); %3rd dimension is FLIM histogram
    plot_FLIM_int(k) = sum(FLIMgroups(k,:));

    intseg = int_masks(:,:,k).* imrotate(image_stack,-mask_angle_rot,rot_method,'crop');
    image_stack_temp = intseg;
    image_stack_temp(image_stack_temp==0)=NaN;
    %int_groups(k) = nansum(image_stack_temp(:));
    int_groups_phocount(k) = nanmean(image_stack_temp(:))...
        *sum(sum(int_masks(:,:,k)))*num_images;
     int_groups(k) = nanmean(image_stack_temp(:));
end

figure(4); clf; plot(int_groups_phocount, plot_FLIM_int, 'bo');
title('Photon counts match between intensity image and FLIMage?');
xlabel('intensity image photons per mask');
ylabel('FLIMage photons per mask'); drawnow;

%%
first_image = int_image0{1};
middle_image = int_image0{round(end/2)};
last_image = int_image0{end};
input_params.mask_type = mask_type;
input_params.scan_mag = scan_mag;
input_params.pixel_bin_width = pixel_bin_width;
input_params.rot_method = rot_method;
reg_seg_plots.first_image = first_image;
reg_seg_plots.middle_image = middle_image;
reg_seg_plots.last_image = last_image;
reg_seg_plots.num_images = num_images;
reg_seg_plots.image_stack = image_stack;
reg_seg_plots.offset_corrected_image_stack = imrotate(image_stack,...
    -mask_angle_rot,rot_method,'crop');
reg_seg_plots.int_masks = int_masks;
reg_seg_plots.mask_distance = mask_distance;
reg_seg_plots.input_params = input_params;
beep;
end




