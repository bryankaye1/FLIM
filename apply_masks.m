function [int_groups_phocount,plot_FLIM_int,int_groups,FLIMgroups] = ...
    apply_masks(int_masks,FLIM_stack,image_stack,mask_angle_rot,rot_method,num_images)

[~,~,num_bins] = size(FLIM_stack);

for k = 1:size(int_masks,3)
    FLIM_masks = repmat(int_masks(:,:,k),1,1,num_bins);
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

b=1; %test tmax-tmin+1 = num_bins
end