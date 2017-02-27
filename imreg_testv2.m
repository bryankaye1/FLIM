%This file is used to test image segmentation
%Written by Olivia Stiehl and Bryan Kaye
clear;
sim_im = 1;
num_images = 2;%length(dataname_cell);

if ~sim_im
    FOV_pth = 'na';
    FOV_fn = 'none';
    tsm = 0; reach = 0;
    tmini=1; tmaxi = 4096;
    data_path = '/Users/bryankaye/Documents/MATLAB/data/2017-01-26/';
    data_name ='DA_ROUND2_spindle2';
    [tsm,dataname_cell] = find_filenames(data_path,data_name,tsm,1);
    %num_images = length(dataname_cell);
    int_image0 = cell(1,num_images);
    for k = 1:num_images
        [int_image0{k}] = spc2image(tmini,tmaxi,dataname_cell{k},data_path,FOV_fn,FOV_pth,tsm,reach);
    end
else
    data_path = '/Users/bryankaye/Documents/MATLAB/data/sim_images/';
    base_name = 'sim_image';
%    int_image0{1} = imread([data_path,'sim_image1.tif']);
%    int_image0{2} = imread([data_path,'sim_image2.tif']);
     namelist_str = dir([data_path,base_name,'*']);
     for dc_ind = 1:length(namelist_str)
        dataname_cell{dc_ind} = namelist_str(dc_ind).name;
        int_image0{dc_ind} = imread([data_path,dataname_cell{dc_ind}]);
     end
     num_images = length(namelist_str);
end
%% calculate rotation and translation of each frame
%make mask
%find threshold and apply it
for k = 1:num_images
    mean_im0(k) = mean2(int_image0{k});
    std_im0(k) = std2(int_image0{k});
    thresh(k) = mean_im0(k)+0.5*std_im0(k);
    thr_mask{k} = int_image0{k}>thresh(k);
    
    %apply erosion and then dilation
    mi0_pre{k} = imfill(thr_mask{k},'holes');
    se_erode = strel('disk',4,0);
    er_mask{k} = imerode(mi0_pre{k},se_erode);
    se_dialate = strel('disk',4,0);
    mi0{k} = imdilate(er_mask{k},se_dialate);
    %mi0{k}=medfilt2(mi0_pre{k});
    
    s = regionprops(mi0{k},'centroid','MajorAxisLength','MinorAxisLength');
    centroid(1,k) = s.Centroid(1)-65;
    centroid(2,k) = s.Centroid(2)-65;
    mi{k} = imtranslate(mi0{k},[-centroid(1,k),-centroid(2,k)]);
    int_image{k} = im2double(imtranslate(int_image0{k},[-centroid(1,k),-centroid(2,k)]));
    i_bracket{k} = mean(int_image{k}(mi{k}));
end
figure(1); clf; imshow(int_image{1},'InitialMagnification','fit');
figure(2); clf; imshow(int_image{2},'InitialMagnification','fit');
%%
thetamin = -8; thetamax = 8; thetastep = .5;
theta = thetamin:thetastep:thetamax;
ctheta_interpx = thetamin:thetastep/5:thetamax; 

for k = 1:num_images
    di{k} = mi{k} .* ((int_image{k} - i_bracket{k}) ./ i_bracket{k});
end

xmin_pix=1;
xmax_pix=128;
ymin_pix=1;
ymax_pix=128;

for k = 1:num_images-1
    theta_ind = 0;
    for theta_var = theta
        theta_ind = theta_ind + 1;
        term1 = fft2(di{k} ((xmin_pix:xmax_pix),(ymin_pix:ymax_pix)));
        term2 = conj(fft2(imrotate(di{k+1} ((xmin_pix:xmax_pix),...
            (ymin_pix:ymax_pix)),theta_var,'bicubic','crop')));
        term3 = fft2(mi{k} ((xmin_pix:xmax_pix),(ymin_pix:ymax_pix)));
        term4 = conj(fft2(imrotate(mi{k+1} ((xmin_pix:xmax_pix),...
            (ymin_pix:ymax_pix)),theta_var,'bicubic','crop')));
        jansC{k,theta_ind} = ifft2(term1.*term2) ./ ifft2(term3 .* term4);
        test{k,theta_ind}(1:10,1:10) = jansC{k,theta_ind}(1:10,1:10);
        test{k,theta_ind}(11:20,1:10) = jansC{k,theta_ind}(1:10,119:128);
        test{k,theta_ind}(1:10,11:20) = jansC{k,theta_ind}(119:128,1:10);
        test{k,theta_ind}(11:20,11:20) = jansC{k,theta_ind}(119:128,119:128);
    end
end

for k = 1:num_images-1
    for theta_ind = 1:length(theta)
        ctheta(k,theta_ind)=max(test{k,theta_ind}(:));
    end
    ctheta_smooth(k,:)=smooth(ctheta(k,:),9);
    [~,maxind_theta(k)] = max(ctheta_smooth(k,:));
    theta_max_smooth(k)=theta(maxind_theta(k));
       
    ctheta_interpy(k,:) = spline(theta,ctheta_smooth(k,:),ctheta_interpx);
    [~,maxind_theta(k)] = max(ctheta_interpy(k,:));
    theta_max(k)=ctheta_interpx(maxind_theta(k));     
end

%%%%%%%%%%%%% check smoothing \& max position
figure(3); clf;
plot(theta,ctheta,'bo');
hold on;
plot(theta,ctheta_smooth,'r');
plot(ctheta_interpx,ctheta_interpy,'g');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k = 1:num_images-1
    term1 = fft2(di{k} ((xmin_pix:xmax_pix),(ymin_pix:ymax_pix)));
    term2 = conj(fft2(imrotate(di{k+1} ((xmin_pix:xmax_pix),...
        (ymin_pix:ymax_pix)),theta_max(k),'bicubic','crop')));
    term3 = fft2(mi{k} ((xmin_pix:xmax_pix),(ymin_pix:ymax_pix)));
    term4 = conj(fft2(imrotate(mi{k+1} ((xmin_pix:xmax_pix),...
        (ymin_pix:ymax_pix)),theta_max(k),'bicubic','crop')));
    jansC_theta_max{k} = ifft2(term1.*term2) ./ ifft2(term3 .* term4);
    pp{k} = jansC_theta_max{k}(1:10,1:10);
    np{k} = jansC_theta_max{k}(1:10,119:128);
    pn{k} = jansC_theta_max{k}(119:128,1:10);
    nn{k} = jansC_theta_max{k}(119:128,119:128);
end

%%%%%%%%%%%%%%%%%%%%%%%%% FIND LX AND LY

rel_transx(1) = round(-centroid(1,1));
rel_transy(1) = round(-centroid(2,1));
total_transx(1) = sum(rel_transx);
total_transy(1) = sum(rel_transy);

for k = 1:num_images-1
    
    if max(pp{k}(:))> max(pn{k}(:)) && max(pp{k}(:)) > max(np{k}(:)) && max(pp{k}(:)) > max(nn{k}(:))
        [~,maxind_n] = max(pp{k}(:));
        [ly(k),lx(k)] = ind2sub(size(pp{k}),maxind_n);
        mask_transx(k) = lx(k)-1;
        mask_transy(k) = ly(k)-1;
    elseif max(pn{k}(:))> max(pp{k}(:)) && max(pn{k}(:)) > max(np{k}(:)) && max(pn{k}(:)) > max(nn{k}(:))
        [~,maxind_n] = max(pn{k}(:));
        [ly(k),lx(k)] = ind2sub(size(pn{k}),maxind_n);
        mask_transx(k) = lx(k)-1;
        mask_transy(k) = -(10-ly(k)+1);
    elseif max(np{k}(:))> max(pp{k}(:)) && max(np{k}(:)) > max(pn{k}(:)) && max(np{k}(:)) > max(nn{k}(:))
        [~,maxind_n] = max(np{k}(:));
        [ly(k),lx(k)] = ind2sub(size(pn{k}),maxind_n);
        mask_transx(k) = -(10-lx(k)+1);
        mask_transy(k) = ly(k)-1;
    elseif max(nn{k}(:))> max(pp{k}(:)) && max(nn{k}(:)) > max(np{k}(:)) && max(nn{k}(:)) > max(pn{k}(:))
        [~,maxind_n] = max(nn{k}(:));
        [ly(k),lx(k)] = ind2sub(size(nn{k}),maxind_n);
        mask_transx(k) = -(10-lx(k)+1);
        mask_transy(k) = -(10-ly(k)+1);
    end
    
    rel_transx(k+1)= round( mask_transx(k) - ( centroid(1,k+1)-centroid(1,k) ) );
    rel_transy(k+1)= round( mask_transy(k) - ( centroid(2,k+1)-centroid(2,k) ) );
    total_transx(k+1) = sum(rel_transx);
    total_transy(k+1) = sum(rel_transy);
    fprintf('xtrans: %2.1f   ytrans: %2.1f\n', total_transx(k),total_transy(k)); 
    fprintf('theta_max %2.1f \n', theta_max(k));
end
    fprintf('xtrans: %2.1f   ytrans: %2.1f\n', total_transx(k+1),total_transy(k+1)); 

%  [~,maxind_n] = max(test_n{k}(:));
%  [ly(k),lx(k)] = ind2sub(size(test_n{k}),maxind_n);
%     if lx< xlen/2
%         mask_transx(k) = -( lx(k)-1);%added "round()" so that when you get 0, it isn't displayed as -0
%     else
%         mask_transx(k) = xlen-lx(k)+1;
%     end
%     
%     if ly< ylen/2
%         mask_transy(k) = -( ly(k)-1 );
%     else
%         mask_transy(k) = ylen-ly(k)+1;
%     end

%load in FLIM, placing pixels at appropriate place