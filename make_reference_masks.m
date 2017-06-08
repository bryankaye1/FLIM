%make reference spindle masks for NOSPINDLE measurements

%2X = 30638;
%8X = 30635;
%12X = 30640;
i = 30640;

load(['/Users/bryankaye/Documents/MATLAB/data/matin/matin',...
    num2str(i),'.mat']);

mi0 = imfill(seg_results.int_masks(:,:,1),'holes');
int_masks = seg_results.int_masks;
mask_distance = seg_results.mask_distance;
scan_mag = seg_results.input_params.scan_mag;
mask_filename = [input(1,1,1).pth_data,input(1,1,1).dataname];
mask_matin_number = i;


if ~scan_mag==12.8
save_name = ['default_masks_',num2str(scan_mag),'X'];
else
save_name = 'default_masks_12X';
end

save(save_name,'mi0','int_masks','mask_distance','scan_mag',...
    'mask_filename','mask_matin_number');