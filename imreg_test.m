%This file is used to test image segmentation
%Written by Olivia Stiehl and Bryan Kaye 
clear;
FOV_pth = 'na';
FOV_fn = 'none';
tsm = 0; reach = 0;
tmini=1; tmaxi = 4096;

data_pth = '/Users/bryankaye/Documents/MATLAB/data/2017-01-26/';
data_name ='DA_ROUND2_spindle2';
[tsm,dataname_cell] = find_filenames(data_pth,data_name,tsm,1);

%Load in intensity files
tic
num_images = 2;%length(dataname_cell);
int_image0 = cell(1,num_images);
for k = 1:num_images
    [int_image0{k}] = spc2image(tmini,tmaxi,dataname_cell{k},data_pth,FOV_fn,FOV_pth,tsm,reach);
end
toc
%% calculate rotation and translation of each frame
%make mask
%find threshold and apply it
for k = 1:num_images
mean_im0(k) = mean2(int_image0{k});
std_im0(k) = std2(int_image0{k});
thresh(k) = mean_im0(k)+std_im0(k);
%thr_image{k} = int_image{k}(int_image{k}>thresh);
thr_mask{k} = int_image0{k}>thresh(k);
%apply erosion and then dilation
se_erode = strel('disk',2,0);
er_mask = imerode(thr_mask{k},se_erode);
se_dialate = strel('disk',6,0);
mi0{k} = imdilate(er_mask,se_dialate);

s = regionprops(mi0{k}, 'centroid');
centroid(1,k) = 65-s.Centroid(1);
centroid(2,k) = 65-s.Centroid(2);
mi{k} = imtranslate(mi0{k},[centroid(1,k),centroid(2,k)]);
int_image{k} = imtranslate(int_image0{k},[centroid(1,k),centroid(2,k)]);
i_bracket{k} = mean(int_image{k}(mi{k}));
end

figure(1); imshow(mi{1},'InitialMagnification', 'fit');
figure(2); imagesc(int_image{1});
figure(3); imshow(mi{2},'InitialMagnification', 'fit');
figure(4); imagesc(int_image{2});
%%
theta = -5:1:5;
for k = 1:num_images
    di{k} = mi{k} .* ((int_image{k} - i_bracket{k}) ./ i_bracket{k});  
end



for k = 1:num_images-1
theta_ind = 0;
    for theta_var = theta
        theta_ind = theta_ind +1;
        %rot_im{k} = imrotate(di{k},1);
        term1 = fft2(di{k});
        term2 = conj(fft2(imrotate(di{k+1},theta_var,'bicubic','crop')));
        term3 = fft2(mi{k});
        term4 = conj(fft2(imrotate(mi{k+1},theta_var,'bicubic','crop')));
        %jansC{k,theta_ind} = ifft2(term1.*term2) ./ ifft2(term3 .* term4);
        jansC{k,theta_ind} = ifft2(term1.*term2);
    end    
end
theta_indmax = theta_ind;
test = jansC;
for k = 1:num_images-1
    for theta_ind = 1:theta_indmax
      test{k,theta_ind}(test{k,theta_ind}==inf)=-inf;
      [~,maxind] = max(test{k,theta_ind}(:));
      [lx(k,theta_ind),ly(k,theta_ind)] = ind2sub(size(test{k,theta_ind}),maxind);
      tmax(k,theta_ind)=max(max(test{k,theta_ind}));
    end    
end


figure(7);
plot(theta,tmax(1,:));
%rot_im = imrotate(di



%transform frame rotations/translation to each pixels translation




%load in FLIM, placing pixels at appropriate place





% for theta_ind = 1:length(test)
% test{theta_ind}(test{theta_ind}==inf)=0;
% [~,maxind] = max(test{theta_ind}(:));
% [lx(theta_ind),ly(theta_ind)] = ind2sub(size(test{theta_ind}),maxind);
% tmax(theta_ind)=max(max(test{theta_ind}));
% end



