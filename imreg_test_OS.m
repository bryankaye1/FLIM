%This file is used to test image segmentation
%Written by Olivia Stiehl and Bryan Kaye 
clear;
close all;
%FOV_pth = 'na';
%FOV_fn = 'none';
%tsm = 0; reach = 0;
%tmini=1; tmaxi = 4096;

%data_pth = 'C:\Users\BIGSAS\Desktop\FLIM\FLIM-master\';
%data_name ='DA_ROUND2_spindle2';
%[tsm,dataname_cell] = find_filenames(data_pth,data_name,tsm,1);
data_pth = 'C:\Users\BIGSAS\Desktop\FLIM\';
data_name ='initial_image.tif';
data_name_2 ='final_image.tif';
%Load in intensity files
%tic
num_images = 2;%length(dataname_cell);
%int_image0 = cell(1,num_images);
%for k = 1:num_images
%    [int_image0{k}] = spc2image(tmini,tmaxi,dataname_cell{k},data_pth,FOV_fn,FOV_pth,tsm,reach);
%end
%toc

int_image0{1}=imread(data_name);
int_image0{2}=imread(data_name_2);

%% calculate rotation and translation of each frame
%make mask
%find threshold and apply it
for k = 1:num_images
mean_im0(k) = mean2(int_image0{k});
std_im0(k) = std2(int_image0{k});
thresh(k) = mean_im0(k)+0.5*std_im0(k);
%thr_image{k} = int_image{k}(int_image{k}>thresh);
thr_mask{k} = int_image0{k}>thresh(k);
%apply erosion and then dilation
%se_erode = strel('disk',1,0);
%er_mask{k} = imerode(thr_mask{k},se_erode);
%se_dialate = strel('disk',4,0);
%mi0{k} = imdilate(er_mask,se_dialate);
%mi0{k} = imfill(er_mask{k},'holes');
mi0_pre{k} = imfill(thr_mask{k},'holes');
se_erode = strel('disk',4,0);
er_mask{k} = imerode(mi0_pre{k},se_erode);
se_dialate = strel('disk',4,0);
mi0{k} = imdilate(er_mask{k},se_dialate)
%mi0{k}=medfilt2(mi0_pre{k});
%mi0{k}=mi0_pre{k};


figure(1); clf; imagesc(int_image0{1});
axis equal
figure(2); clf; imshow(mi0{1});


%figure('name','processed im 2'); imshow(mi0{2},'InitialMagnification', 'fit');
%figure('name','image int im 2'); imagesc(int_image0{2});
%figure('name','thr_mask1'); imagesc(thr_mask{1});
%figure('name','thr_mask2'); imagesc(thr_mask{2});
%figure('name','2: filled_1'); clf; imshow(mi0_pre{1},'InitialMagnification','fit');
%figure('name','filled_2'); imagesc(mi0_pre{2});


s = regionprops(mi0{k},'centroid','MajorAxisLength','MinorAxisLength');
centroid(1,k) = 65-s.Centroid(1);
centroid(2,k) = 65-s.Centroid(2);
centroid(1,1)=65-s.Centroid(1)+3;
centroid(1,2)=65-s.Centroid(2)+3;
mi{k} = imtranslate(mi0{k},[centroid(1,k),centroid(2,k)]);
int_image{k} = im2double(imtranslate(int_image0{k},[centroid(1,k),centroid(2,k)]));
i_bracket{k} = mean(int_image{k}(mi{k}));
end


%%
theta = 0%-20:0.05:20;
for k = 1:num_images
    di{k} = mi{k} .* ((int_image{k} - i_bracket{k}) ./ i_bracket{k});  
end

%min_pix=64-round(s.MinorAxisLength/2);
%max_pix=64+round(s.MinorAxisLength/2);
min_pix=64-round(s.MinorAxisLength/4);
max_pix=64+round(s.MinorAxisLength/4);

for k = 1:num_images-1
theta_ind = 0;
    for theta_var = theta
        theta_ind = theta_ind +1;
        %rot_im{k} = imrotate(di{k},1);
        term1 = fft2(di{k} ((min_pix:max_pix),(min_pix:max_pix)));
        term2 = conj(fft2(imrotate(di{k+1} ((min_pix:max_pix),(min_pix:max_pix)),theta_var,'bicubic','crop')));
        term3 = fft2(mi{k} ((min_pix:max_pix),(min_pix:max_pix)));
        term4 = conj(fft2(imrotate(mi{k+1} ((min_pix:max_pix),(min_pix:max_pix)),theta_var,'bicubic','crop')));
        jansC{k,theta_ind} = ifft2(term1.*term2) ./ ifft2(term3 .* term4);
        %jansC{k,theta_ind} = ifft2(term1.*term2);
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
    tmax_smooth(k,:)=smooth(tmax(k,:),9);
    
    [~,maxind_theta(k)] = max(tmax_smooth(k,:));
    theta_max(k)=theta(maxind_theta(k));
    
end




%%%%%%%%%%%%% check smoothing \& max position
figure(3);
plot(theta,tmax(1,:));
hold on;
plot(theta,tmax_smooth,'r');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
for k = 1:num_images-1

        term1 = fft2(di{k} ((min_pix:max_pix),(min_pix:max_pix)));
        term2 = conj(fft2(imrotate(di{k+1} ((min_pix:max_pix),(min_pix:max_pix)),theta_max(k),'bicubic','crop')));
        term3 = fft2(mi{k} ((min_pix:max_pix),(min_pix:max_pix)));
        term4 = conj(fft2(imrotate(mi{k+1} ((min_pix:max_pix),(min_pix:max_pix)),theta_max(k),'bicubic','crop')));
        jansC_theta_max{k} = ifft2(term1.*term2) ./ ifft2(term3 .* term4);
         
end

%%%%%%%%%%%%%%%%%%%%%%%%% FIND LX AND LY

test_n = jansC_theta_max;
for k = 1:num_images-1
   
      test_n{k}(test_n{k}==inf)=-inf;
      [~,maxind_n] = max(test_n{k}(:));
      [lx_n(k),ly_n(k)] = ind2sub(size(test_n{k}),maxind_n);
      %tmax(k)=max(max(test{k,theta_ind}));
   
    %tmax_smooth(k,:)=smooth(tmax(k,:),9);
    
    %[~,maxind_theta(k)] = max(tmax_smooth(k,:));
    %theta_max(k)=theta(maxind_theta(k));
    %translation_x_final(k)=round(lx_n(k)+centroid(1,k));
    %translation_y_final(k)=round(ly_n(k)+centroid(2,k));
    
end



%load in FLIM, placing pixels at appropriate place

