%Script to take a image, translate and rotate it, and then add gaussian
%noise.
rotation_degs = {0,1,2,3,2,1,0};
translation_pix = {[0,0],[1,1],[2,2],[3,3],[4,4],[5,5],[6,6]};
noise_mean = 0;
noise_var = 0;

%% Load in image
FOV_pth = 'na';
FOV_fn = 'none';
tsm = 0; reach = 0;
tmini=1; tmaxi = 4096;

data_pth = '/Users/bryankaye/Documents/MATLAB/data/2017-01-26/';
data_name ='DA_ROUND2_spindle2_c01.sdt';
[raw_image] = spc2image(tmini,tmaxi,data_name,data_pth,FOV_fn,FOV_pth,tsm,reach);
initial_image = mat2gray(raw_image,[0,max(max(raw_image))]);

%% Translate, Rotate Image, Add noise, and save image
for i = 1:length(rotation_degs)
    rot_im{i} = imrotate(initial_image,rotation_degs{i},'crop');
    rot_tr_im{i} = imtranslate(rot_im{i},translation_pix{i});
    sim_image{i} = imnoise(rot_tr_im{i},'gaussian',noise_mean,noise_var);
    imwrite(sim_image{i},['/Users/bryankaye/Documents/MATLAB/data/sim_images/sim_image',...
        num2str(i),'.tif']);
end
figure(1); clf; 
subplot(1,3,1); imshow(initial_image,'InitialMagnification','fit');
subplot(1,3,2); imshow(sim_image{1},'InitialMagnification','fit');
subplot(1,3,3); imshow(sim_image{end},'InitialMagnification','fit');