%This file is used to test image segmentation
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
int_image = cell(1,num_images);
for k = 1:num_images
    [int_image{k}] = spc2image(tmini,tmaxi,dataname_cell{k},data_pth,FOV_fn,FOV_pth,tsm,reach);
end
toc
%% calculate rotation and translation of each frame
%make mask
%find threshold and apply it
for k = 1:num_images
mean_im(k) = mean2(int_image{k});
std_im(k) = std2(int_image{k});
thresh(k) = mean_im(k)+std_im(k);
%thr_image{k} = int_image{k}(int_image{k}>thresh);
thr_mask{k} = int_image{k}>thresh(k);
%apply erosion and then dilation
se_erode = strel('disk',2,0);
er_mask = imerode(thr_mask{k},se_erode);
se_dialate = strel('disk',6,0);
mi{k} = imdilate(er_mask,se_dialate);
end

figure(1); imshow(mi{1},'InitialMagnification', 'fit');
figure(2); imagesc(int_image{1});
figure(1); imshow(mi{2},'InitialMagnification', 'fit');
figure(2); imagesc(int_image{2});
%%
theta = -10:1:-1;
for k = 1:num_images
    di{k} = thr_mask{k} .* ((int_image{k} - mean_im(k)) ./ mean_im(k));  
end

tic
ind = 0;
for theta_var = theta
    ind = ind +1;
    for k = 1:num_images-1
        %rot_im{k} = imrotate(di{k},1);
        term1 = fft2(di{k});
        term2 = conj(fft2(imrotate(di{k+1},theta_var,'crop')));
        term3 = fft2(mi{k});
        term4 = conj(fft2(imrotate(mi{k+1},theta_var,'crop')));
        jansC{ind} = ifft2(term1.*term2) ./ ifft2(term3 .* term4);
    end
    
end
toc


test=jansC;
for ind = 1:length(test)
test{ind}(test{ind}==inf)=0;
[~,maxind] = max(test{ind}(:));
[lx(ind),ly(ind)] = ind2sub(size(test{ind}),maxind);
tmax(ind)=max(max(test{ind}));
end

figure(7);
plot(theta,tmax);
%rot_im = imrotate(di



%transform frame rotations/translation to each pixels translation




%load in FLIM, placing pixels at appropriate place




