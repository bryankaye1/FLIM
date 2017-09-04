%This file is used to test image segmentation
%Written by Olivia Stiehl and Bryan Kaye

function [FLIM_region,ni_region,ngr, segment_results] = make_regions(data_path,data_name,...
    tmini,tmaxi,FOV_fn,FOV_pth,mask_type,scan_mag,pixel_bin_width,angle_dep)

%[FLIM_region,ni_region,ngr, segment_resutls]
%[FLIMgroups,int_groups,reg_seg_plots]

sim_im = 1;
rot_method = 'nearest'; %just did bilinear
if ~sim_im
    %FOV_pth = 'na'; FOV_fn = 'none'; mask_type = 'ellipsoid'; scan_mag = 8;
    %tmini=1; tmaxi = 4096; data_path = '/Users/bryankaye/Documents/MATLAB/data/2017-01-26/';
    %data_name ='Donor_r1'; 
    [tsm,dataname_cell] = find_filenames(data_path,[data_name,'_c'],0,1);
    num_images = length(dataname_cell);

    for k = 1:num_images     
        [int_image0{k},FLIMages{k}] = spc2FLIMdata(tmini,tmaxi,...
            dataname_cell{k},data_path,FOV_fn,FOV_pth,tsm,0);
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


if ~contains(data_name,'NOSPINDLE')
    %Make image_stack by registering images
    [image_stack,registration_vectors] = register_images(int_image0,...
        scan_mag,rot_method);
    
    %Make masks from registered images
    if strcmp(mask_type,'ellipsoid')
        [int_masks,mask_distance,mask_angle_rot] = ellips_seg(image_stack,...
            pixel_length,scan_mag,rot_method);
    elseif strcmp(mask_type, 'edge_distance')
        [int_masks,mask_distance,mask_angle_rot,mi0] = dist_seg(int_image0,...
            image_stack,pixel_length,pixel_bin_width,scan_mag);
    else
        input('ERROR: NO MASK TYPE SELECTED');
    end
    
    %Make FLIM_stack by registering FLIM data.
    [FLIM_stack] = transform_FLIMage(FLIMages,registration_vectors,...
        mask_angle_rot,rot_method);
else
    
    %make image_stack
    image_stack = 0;
    for k=1:num_images
        image_stack = int_image0{k}+image_stack;
    end
    
    %load in int_masks and mask_distance of a real spindle image
    if ~scan_mag==12.8
    load(['default_masks_',num2str(scan_mag),'X']);
    else
    load('default_masks_12X');
    end
    mask_angle_rot = 0;
    
    %make FLIM_stack
    FLIM_stack = 0;
    for k = 1:length(FLIMages)
        FLIM_stack = FLIM_stack + FLIMages{k};
    end
    %find better names for vars
    reg_seg_plots.mask_filename = mask_filename;
    reg_seg_plots.mask_matin_number = mask_matin_number;
end


%%
% [int_groups_phocount,plot_FLIM_int,int_groups,FLIMgroups] = ...
%     apply_masks(int_masks,FLIM_stack,image_stack,mask_angle_rot,...
%     rot_method,num_images);

dist_in = -3;
dist_out = 2;

%combine masks into three regions

%make masks for each region
masks_index_in = find(mask_distance<dist_in);

masks_t1 = find(mask_distance>=dist_in);
masks_t2 = find(mask_distance<=dist_out);
masks_index_boundary = intersect(masks_t1,masks_t2);
masks_index_out = find(mask_distance>dist_out);

masks_in = sum(int_masks(:,:,masks_index_in),3);
masks_boundary = sum(int_masks(:,:,masks_index_boundary),3);
masks_out = sum(int_masks(:,:,masks_index_out),3);
masks = {masks_in,masks_boundary,masks_out};

for i = 1:length(masks)
    
    FLIM_rs = reshape(FLIM_stack,[],size(FLIM_stack,3),1);
    masks_rs = reshape(FLIM_stack,[],size(FLIM_stack,3),1);
    
    
    FLIM_masks = repmat(masks{i},1,1,num_bins);
    FLIMseg = FLIM_masks .* FLIM_stack;
    FLIMseg_rs = reshape(A,[],size(A,3),1);
    %FLIMgroups(k,:) = squeeze(sum(sum(FLIMseg))); %3rd dimension is FLIM histogram
    
end
figure(4); clf; plot(int_groups_phocount, plot_FLIM_int, 'bo');
title('Photon counts match between intensity image and FLIMage?');
xlabel('intensity image photons per mask');
ylabel('FLIMage photons per mask'); drawnow;


if angle_dep
    mi0_props = regionprops(mi0,'centroid','MajorAxisLength',...
        'MinorAxisLength','Orientation');
    spindle_angle=mi0_props.Orientation;
    
    zero128 = zeros(128,128);
    one128 = ones(128,128);
    quad{1} = [one128,zero128;zero128,zero128];
    quad{2} = [zero128,one128;zero128,zero128];
    quad{3} = [zero128,zero128;one128,zero128];
    quad{4} = [zero128,zero128;zero128,one128];
    
    for i = 1:length(quad)
        rot_quad{i} = imrotate(quad{i},45+spindle_angle,'nearest','crop');
        rs_rot_quad{i} = rot_quad{i}(65:192,65:192);
        int_masks_quad = int_masks.*rs_rot_quad{i};
        [int_groups_phocount,plot_FLIM_int,qint_groups,qFLIMgroups] = ...
            apply_masks(int_masks_quad,FLIM_stack,image_stack,...
            mask_angle_rot,rot_method,num_images);     

        int_groups = [int_groups,qint_groups];
        FLIMgroups = [FLIMgroups;qFLIMgroups]; %Concatanate in a way that works for FLIM histrogra,s
        seg_ind{i} = 1+(i-1)*length(qint_groups):i*length(qint_groups);
        quad_mask{i} = int_masks_quad;
    end  
  reg_seg_plots.div_groups = seg_ind; 
  reg_seg_plots.div_masks = quad_mask; 
end

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



% for k = 1:size(int_masks,3)
%     FLIM_masks = repmat(int_masks(:,:,k),1,1,tmaxi-tmini+1);
%     FLIMseg = FLIM_masks .* FLIM_stack;
%     FLIMgroups(k,:) = squeeze(sum(sum(FLIMseg))); %3rd dimension is FLIM histogram
%     plot_FLIM_int(k) = sum(FLIMgroups(k,:));
% 
%     intseg = int_masks(:,:,k).* imrotate(image_stack,-mask_angle_rot,rot_method,'crop');
%     image_stack_temp = intseg;
%     image_stack_temp(image_stack_temp==0)=NaN;
%     %int_groups(k) = nansum(image_stack_temp(:));
%     int_groups_phocount(k) = nanmean(image_stack_temp(:))...
%         *sum(sum(int_masks(:,:,k)))*num_images;
%      int_groups(k) = nanmean(image_stack_temp(:));
% end

