%This file is used to test image segmentation
%Written by Olivia Stiehl and Bryan Kaye

function [FLIM_stack,image_stack,reg_seg_plots] = register_FLIMages(data_path,data_name,...
    tmini,tmaxi,FOV_fn,FOV_pth,scan_mag)


sim_im = 0;
rot_method = 'bicubic';
if ~sim_im
    %FOV_pth = 'na'; FOV_fn = 'none'; mask_type = 'ellipsoid'; scan_mag = 8; tmini=1; tmaxi = 4096; data_path = '/Users/bryankaye/Documents/MATLAB/data/2017-01-26/'; data_name ='Donor_r1'; 
    tsm = 0; reach = 0;
    [tsm,dataname_cell] = find_filenames(data_path,data_name,tsm,1);
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

[FLIM_stack] = transform_FLIMage(FLIMages,registration_vectors,...
    0,rot_method);
%%
first_image = int_image0{1};
middle_image = int_image0{round(end/2)};
last_image = int_image0{end};

reg_seg_plots.first_image = first_image;
reg_seg_plots.middle_image = middle_image;
reg_seg_plots.last_image = last_image;
reg_seg_plots.num_images = num_images;
reg_seg_plots.image_stack = image_stack;
reg_seg_plots.offset_corrected_image_stack = imrotate(image_stack,...
    0,rot_method,'crop');

end




