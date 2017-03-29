function [image_stack,registration_vectors] = register_images(int_image0,...
    scan_mag,rot_method)
%This function takes in a cell array of 128x128 images and registers them
%accoring to the algorithm described by Jan Brugues.
%It translates images by up to mcd (max convolution distance)
%in any dimension and rotates by up to 180 degrees.
%mcd is max convolutoin distance
%Code written by Olivia Stiehl and Bryan Kaye 2017.
num_images = length(int_image0);
% show_raw = 1;
% if show_raw
%     for k = 1:num_images
%         figure(4+k); clf;
%         imshow(mat2gray(int_image0{k}));
%     end
% end
% drawnow;

for k = 1:num_images
%     %mean_im0(k) = mean2(int_image0{k});
%     %std_im0(k) = std2(int_image0{k});
%     %thresh(k) = mean_im0(k)+(4/scan_mag)*std_im0(k);
%     thresh(k) = mean2(int_image0{k}(1:3,1:3))*1.5;
%     thr_mask{k} = int_image0{k}>thresh(k);
%     
%     %apply erosion and then dilation
%     mi0_pre{k} = imfill(thr_mask{k},'holes');
%     se_erode = strel('disk',scan_mag,0);
%     er_mask{k} = imerode(mi0_pre{k},se_erode);
%     se_dialate = strel('disk',scan_mag,0);
%     mi0{k} = imdilate(er_mask{k},se_dialate);
%     mi0{k} = imfill(mi0{k},'holes');
    %mi0{k}=medfilt2(mi0_pre{k});

    [mi0{k}] = make_masks(int_image0{k},floor(scan_mag),1.5,1:3);
    s = regionprops(mi0{k},'centroid','MajorAxisLength','MinorAxisLength');
    
    %IF MULTIPLE REGOINS, REMOVE REGIONS SMALLER THAN THRESHOLD AREA
    %IF ZERO REGIONS, REDUCE THRESHOLD
    if length(s)~=1
        [mi0{k}] = make_masks(int_image0{k},floor(scan_mag),1.25,1:3);
        s = regionprops(mi0{k},'centroid','MajorAxisLength','MinorAxisLength');
    end
    
    if length(s)~=1
        [mi0{k}] = make_masks(imgaussfilt(int_image0{k},2),floor(scan_mag),1.5,1:3);
        if scan_mag==2
            mi0{k}(90:128,:)=0;
            mi0{k}(1:30,:)=0;
            mi0{k}(:,90:128)=0;
            mi0{k}(:,1:30)=0;
        end
        s = regionprops(mi0{k},'centroid','MajorAxisLength','MinorAxisLength');
    end
    
    if length(s)~=1
        [mi0{k}] = make_masks(imgaussfilt(int_image0{k},2),floor(scan_mag),2,1:3);
        if scan_mag==2
            mi0{k}(90:128,:)=0;
            mi0{k}(1:30,:)=0;
            mi0{k}(:,90:128)=0;
            mi0{k}(:,1:30)=0;
        end
        s = regionprops(mi0{k},'centroid','MajorAxisLength','MinorAxisLength');
end
    
    centroid(1,k) = s.Centroid(1)-65;
    centroid(2,k) = s.Centroid(2)-65;
    mi{k} = imtranslate(mi0{k},[-centroid(1,k),-centroid(2,k)]);
    int_image{k} = im2double(imtranslate(int_image0{k},[-centroid(1,k),-centroid(2,k)]));
    i_bracket{k} = mean(int_image{k}(mi{k}));
end
%figure(1); clf; imshow(int_image{1},'InitialMagnification','fit');
%figure(2); clf; imshow(int_image{2},'InitialMagnification','fit');
%%
thetamin = -8; thetamax = 8; thetastep = .5;
theta = thetamin:thetastep:thetamax;
ctheta_interpx = thetamin:thetastep/5:thetamax;

for k = 1:num_images
    di{k} = mi{k} .* ((int_image{k} - i_bracket{k}) ./ i_bracket{k});
end

figure(1); clf;
subplot(1,2,1);
imshow(mat2gray(di{1}),'InitialMagnification','fit');
title('first normalized masked image');

subplot(1,2,2);
imshow(mat2gray(di{end}),'InitialMagnification','fit');
title('last normalized masked image');

%This was used to resrtict the pixels that are used to find maximum
%correlation. Due to periodic boundary conditions in the discrete fourier
%transform, we want to use entire image in the fourier transforms.
xmin_pix=1; xmax_pix=128; ymin_pix=1; ymax_pix=128;

short_axis_len = 10 ./ ((440./scan_mag) / 128); %minimum spindle short axis length in pixels 
mcd = floor(short_axis_len*.5); %max convolution distance



for k = 1:num_images-1
    theta_ind = 0;
    for theta_var = theta
        theta_ind = theta_ind + 1;
        term1 = fft2(di{k} ((xmin_pix:xmax_pix),(ymin_pix:ymax_pix)));
        term2 = conj(fft2(imrotate(di{k+1} ((xmin_pix:xmax_pix),...
            (ymin_pix:ymax_pix)),theta_var,rot_method,'crop')));
        term3 = fft2(mi{k} ((xmin_pix:xmax_pix),(ymin_pix:ymax_pix)));
        term4 = conj(fft2(imrotate(mi{k+1} ((xmin_pix:xmax_pix),...
            (ymin_pix:ymax_pix)),theta_var,rot_method,'crop')));
        jansC{k,theta_ind} = ifft2(term1.*term2) ./ ifft2(term3 .* term4);       
        test{k,theta_ind}(1:mcd,1:mcd) = jansC{k,theta_ind}(1:mcd,1:mcd);
        test{k,theta_ind}(mcd+1:2*mcd,1:mcd) = jansC{k,theta_ind}(1:mcd,128-mcd+1:128);
        test{k,theta_ind}(1:mcd,mcd+1:2*mcd) = jansC{k,theta_ind}(128-mcd+1:128,1:mcd);
        test{k,theta_ind}(mcd+1:2*mcd,mcd+1:2*mcd) = ...
            jansC{k,theta_ind}(128-mcd+1:128,128-mcd+1:128);
    end
end

for k = 1:num_images-1
    for theta_ind = 1:length(theta)
        ctheta(k,theta_ind)=max(test{k,theta_ind}(:));
    end
    ctheta_smooth(k,:)=smooth(ctheta(k,:),11);
    [~,maxind_theta(k)] = max(ctheta_smooth(k,:));
    theta_max_smooth(k)=theta(maxind_theta(k));
    
    ctheta_interpy(k,:) = spline(theta,ctheta_smooth(k,:),ctheta_interpx);
    [~,maxind_theta(k)] = max(ctheta_interpy(k,:));
    theta_max(k)=ctheta_interpx(maxind_theta(k));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k = 1:num_images-1
    term1 = fft2(di{k} ((xmin_pix:xmax_pix),(ymin_pix:ymax_pix)));
    term2 = conj(fft2(imrotate(di{k+1} ((xmin_pix:xmax_pix),...
        (ymin_pix:ymax_pix)),theta_max(k),rot_method,'crop')));
    term3 = fft2(mi{k} ((xmin_pix:xmax_pix),(ymin_pix:ymax_pix)));
    term4 = conj(fft2(imrotate(mi{k+1} ((xmin_pix:xmax_pix),...
        (ymin_pix:ymax_pix)),theta_max(k),rot_method,'crop')));
    jansC_theta_max{k} = ifft2(term1.*term2) ./ ifft2(term3 .* term4);
    pp{k} = jansC_theta_max{k}(1:mcd,1:mcd);
    np{k} = jansC_theta_max{k}(1:mcd,128-mcd+1:128);
    pn{k} = jansC_theta_max{k}(128-mcd+1:128,1:mcd);
    nn{k} = jansC_theta_max{k}(128-mcd+1:128,128-mcd+1:128);
end

image_rotation(1) = 0;
mask_cumulative_transx(1) = 0;
mask_cumulative_transy(1) = 0;
%% Record transformations vectors
for k = 1:num_images-1
    %This section calculates how much the mask has to move after it's been
    %centered and rotated by theta_max. Also keep tracks of the angle each
    %image has to rotate by
    if max(pp{k}(:))> max(pn{k}(:)) && max(pp{k}(:)) > max(np{k}(:)) && max(pp{k}(:)) > max(nn{k}(:))
        [~,maxind_n] = max(pp{k}(:));
        [ly(k),lx(k)] = ind2sub(size(pp{k}),maxind_n);
        mask_transx(k) = lx(k)-1;
        mask_transy(k) = ly(k)-1;
    elseif max(pn{k}(:))> max(pp{k}(:)) && max(pn{k}(:)) > max(np{k}(:)) && max(pn{k}(:)) > max(nn{k}(:))
        [~,maxind_n] = max(pn{k}(:));
        [ly(k),lx(k)] = ind2sub(size(pn{k}),maxind_n);
        mask_transx(k) = lx(k)-1;
        mask_transy(k) = -(mcd-ly(k)+1);
    elseif max(np{k}(:))> max(pp{k}(:)) && max(np{k}(:)) > max(pn{k}(:)) && max(np{k}(:)) > max(nn{k}(:))
        [~,maxind_n] = max(np{k}(:));
        [ly(k),lx(k)] = ind2sub(size(pn{k}),maxind_n);
        mask_transx(k) = -(mcd-lx(k)+1);
        mask_transy(k) = ly(k)-1;
    elseif max(nn{k}(:))> max(pp{k}(:)) && max(nn{k}(:)) > max(np{k}(:)) && max(nn{k}(:)) > max(pn{k}(:))
        [~,maxind_n] = max(nn{k}(:));
        [ly(k),lx(k)] = ind2sub(size(nn{k}),maxind_n);
        mask_transx(k) = -(mcd-lx(k)+1);
        mask_transy(k) = -(mcd-ly(k)+1);
    end    
    image_rotation(k+1) = image_rotation(k) + theta_max(k);
    mask_cumulative_transx(k+1) = sum(mask_transx);
    mask_cumulative_transy(k+1) = sum(mask_transy);    
    %   fprintf('xtrans: %2.1f   ytrans: %2.1f\n', total_transx(k),total_transy(k));
    %   fprintf('theta_max %2.1f \n', theta_max(k));
end

% if any([mask_cumulative_transx,mask_cumulative_transy])
%     fprintf('mask was moved after rotation\n');
% end

registration_vectors.rotation = image_rotation;
registration_vectors.centroid = centroid;
registration_vectors.last_translationx = mask_cumulative_transx;
registration_vectors.last_translationy = mask_cumulative_transy;

%% Transform Images and then align them

for k = 1:num_images
    translated_image{k} = imtranslate(int_image0{k},[-centroid(1,k),-centroid(2,k)]);
    tran_rot_image{k} = imrotate(translated_image{k},image_rotation(k),rot_method,'crop');
    tran_rot_tran_image{k} = imtranslate(tran_rot_image{k},[mask_cumulative_transx(k),mask_cumulative_transy(k)]);
end

% image_stack = 0;
% for k = 1:num_images
%     image_stack = image_stack + tran_rot_tran_image{k};    
% end

image_stack_3D = zeros(128,128,num_images);
tran_rot_tran_image_NaN = tran_rot_tran_image;
 for k = 1:num_images
     tran_rot_tran_image_NaN{k}(tran_rot_tran_image_NaN{k}==0) = NaN;
     image_stack_3D(:,:,k) = tran_rot_tran_image_NaN{k};    
 end
 
image_stack = nanmean(image_stack_3D,3);
image_stack(isnan(image_stack))=0;


h2 = figure(2); clf;
subplot(2,2,1); imshow(mat2gray(int_image{1}),'InitialMagnification','fit');
title('first intensity image');
subplot(2,2,2); imshow(mat2gray(int_image{round(num_images/2)}),'InitialMagnification','fit');
title('middle intensity image');
subplot(2,2,3); imshow(mat2gray(int_image{end}),'InitialMagnification','fit');
title('last intensity image');
subplot(2,2,4); imshow(mat2gray(image_stack),'InitialMagnification','fit');
title('registered intensity image');
figure(5); clf; subplot(3,ceil(num_images/3),1);

for k = 1:num_images
    subplot(3,ceil(num_images/3),k);
    imshow(mat2gray(int_image{k}),'InitialMagnification','fit');
end
drawnow;

end

