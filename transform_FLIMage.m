
function [FLIM_stack] = transform_FLIMage(FLIMages,registration_vectors,...
    mask_angle_rot,rot_method)

centroid = registration_vectors.centroid;
image_rotation = registration_vectors.rotation;
mask_cumulative_transx = registration_vectors.last_translationx;
mask_cumulative_transy = registration_vectors.last_translationy;

FLIM_stack = 0;

for k = 1:length(FLIMages)
    
    translated_FLIMage = imtranslate(FLIMages{k},[-centroid(1,k),...
        -centroid(2,k)]);
    tran_rot_FLIMage = imrotate(translated_FLIMage,image_rotation(k),...
        rot_method,'crop');
    tran_rot_tran_FLIMage = imtranslate(tran_rot_FLIMage,...
        [mask_cumulative_transx(k),mask_cumulative_transy(k)]); 
    final_FLIMage = imrotate(tran_rot_tran_FLIMage,-mask_angle_rot,...
        rot_method,'crop');
    
    FLIM_stack = FLIM_stack + final_FLIMage;

end



end