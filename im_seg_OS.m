%This file is used to test image segmentation
%Written by Olivia Stiehl and Bryan Kaye 
clear;
%close all;

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


s{k} = regionprops(mi0{k},'centroid','MajorAxisLength','MinorAxisLength','Orientation');
% calculate vectors how center must be shifted to reach center of image
centroid(1,k) = 65-s{k}.Centroid(1);
centroid(2,k) = 65-s{k}.Centroid(2);
mi{k} = imtranslate(mi0{k},[centroid(1,k),centroid(2,k)]);
int_image{k} = im2double(imtranslate(int_image0{k},[centroid(1,k),centroid(2,k)]));
i_bracket{k} = mean(int_image{k}(mi{k}));
end


%%
theta = -20:0.05:20;
for k = 1:num_images
    di{k} = mi{k} .* ((int_image{k} - i_bracket{k}) ./ i_bracket{k});  

end

% ASSUME THAT SPINDLE SIZE DOES NOT CHANGE DRASTICALLY OVER TIME OF
% EXPERIMENT
min_pix=64-round(s{1}.MinorAxisLength/2);
max_pix=64+round(s{1}.MinorAxisLength/2);


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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% ACTUALLY WE ASSUME FIRST TRANSLATION IS PRECISE ENOUGH
%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%for k = 1:num_images-1

 %       term1 = fft2(di{k} ((min_pix:max_pix),(min_pix:max_pix)));
 %       term2 = conj(fft2(imrotate(di{k+1} ((min_pix:max_pix),(min_pix:max_pix)),theta_max(k),'bicubic','crop')));
 %       term3 = fft2(mi{k} ((min_pix:max_pix),(min_pix:max_pix)));
 %       term4 = conj(fft2(imrotate(mi{k+1} ((min_pix:max_pix),(min_pix:max_pix)),theta_max(k),'bicubic','crop')));
 %       jansC_theta_max{k} = ifft2(term1.*term2) ./ ifft2(term3 .* term4);
         
%end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% FIND LX AND LY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%test_n = jansC_theta_max;
%for k = 1:num_images-1
   
%      test_n{k}(test_n{k}==inf)=-inf;
%      [~,maxind_n] = max(test_n{k}(:));
%      [lx_n(k),ly_n(k)] = ind2sub(size(test_n{k}),maxind_n);
     
   
    %theta_max(k)=theta(maxind_theta(k));
    %translation_x_final(k)=round(lx_n(k)+centroid(1,k));
    %translation_y_final(k)=round(ly_n(k)+centroid(2,k));
    
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

angle_offset=s{1}.Orientation;
%%%%%%%%%%%%%%%%%%% in order to check this offset - create overlay first
%%%%%%%%%%%%%%%%%%% image, elliposid, MajorAxis and horizontal image axis;
figure(4)
imagesc(mi0{1});
hold on;

phi_t = linspace(0,2*pi,50);
cosphi_t = cos(phi_t);
sinphi_t = sin(phi_t);


    xbar_t = s{1}.Centroid(1);
    ybar_t = s{1}.Centroid(2);

    a_t = s{1}.MajorAxisLength/2;
    b_t = s{1}.MinorAxisLength/2;

    theta_t = pi*s{1}.Orientation/180;
    R_t = [ cos(theta_t)   sin(theta_t)
         -sin(theta_t)   cos(theta_t)];

    xy_t = [a_t*cosphi_t; b_t*sinphi_t];
    xy_t = R_t*xy_t;

    x_t = xy_t(1,:) + xbar_t;
    y_t = xy_t(2,:) + ybar_t;

    plot(x_t,y_t,'r','LineWidth',2);

hold on

hlen = s{1}.MajorAxisLength/2;
cosOrient_tt = cosd(s{1}.Orientation);
sinOrient_tt = sind(s{2}.Orientation);
xcoords = xbar_t + hlen * [cosOrient_tt -cosOrient_tt];
ycoords = ybar_t + hlen * [-sinOrient_tt sinOrient_tt];
line(xcoords, ycoords);



x_ax=round(s{1}.Centroid-s{1}.MajorAxisLength):1:round(s{1}.Centroid+s{1}.MajorAxisLength);
hold on;
plot(s{1}.Centroid(1),s{1}.Centroid(2),'kx');
hold on;
y_hor(1:length(x_ax))=s{1}.Centroid(2);
plot(x_ax,y_hor,'k');
hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%% VALUES to rotate: theta_max + angle_offset from first image!
int_final{1}=imrotate(int_image{k} ,(-angle_offset),'bicubic','crop');
mi_final{1}=imrotate(mi{k} ,(-angle_offset),'bicubic','crop');

for k=1:num_images-1
int_final{k+1}=imrotate(int_image{k+1} ,(theta_max(k)-angle_offset),'bicubic','crop');
mi_final{k+1}=imrotate(mi{k+1} ,(theta_max(k)-angle_offset),'bicubic','crop');
%figure(6)
%imagesc(mi_final{k+1});
end


for k=1:num_images
figure
imagesc(int_final{k});   
end    


disp('Confirm image alignment by pressing any key !')
pause;
%%%%% HERE: all spindles are centered and aligned along the image axis;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% next: create ring masks



for k=1:num_images
    
mag_fac=1:0.25:3;   

for i=1:length(mag_fac);
    
    %%% only allow for ellipses that fit the image:
    if 65-(mag_fac(i)*s{k}.MajorAxisLength/2) < 1 | 65 + (mag_fac(i)*s{k}.MajorAxisLength/2) > size(int_final{k} (:,2)) | 65-(mag_fac(i)*s{k}.MinorAxisLength/2) < 1 | 65 + (mag_fac(i)*s{k}.MinorAxisLength/2) > size(int_final{k} (:,1))
    mag_fac(i)=NaN;
    else
    %%% create masks of the respective ellipses for each images:
    [X Y]= meshgrid(1:128,1:128); 
    ellipse{k}(:,:,i) = ((X-65)/(mag_fac(i)*s{k}.MajorAxisLength/2)).^2+((Y-65)/(mag_fac(i)*s{k}.MinorAxisLength/2)).^2<=1; 
    end
    
end
    %figure
    %imshow(ellipse{1}(:,:,5));
    %%% create rings:
        for i=1:(size(ellipse{k},3)-1);
        rings{k} (:,:,i) = ellipse{k}(:,:,i+1)-ellipse{k}(:,:,i);
        end    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% ATTEMPT TO VISUALIZE POSITION OF RINGS ON TOP OF THE INTENSITY 
for k=1:num_images
figure
imagesc(int_final{k})
axis equal
hold on;


for i=1:(size(ellipse{k},3)-1);
%phi_t = linspace(0,2*pi,50);
%cosphi_t = cos(phi_t);
%sinphi_t = sin(phi_t);



    %a_t = s{k}.MajorAxisLength/2;
    %b_t = s{1}.MinorAxisLength/2;

    %theta_t = pi*s{1}.Orientation/180;
    R_norm = [ 1   0
            0   1];

    xy_norm = [mag_fac(i)*s{k}.MajorAxisLength/2*cosphi_t; mag_fac(i)*s{k}.MinorAxisLength/2*sinphi_t];
    xy_norm = R_norm*xy_norm;

    x_norm(:,:,i) = xy_norm(1,:) + 65;
    y_norm(:,:,i) = xy_norm(2,:) + 65;

    plot(x_norm(:,:,i),y_norm(:,:,i),'r','LineWidth',2);
    hold on;
end
end


%%%%%%%%%%%%%%%%%%%%%%% OVERLAY THESE RINGS WITH IMAGES;
%%%%%%%%%%% IN CASE OF MEAN - SET "res" TO NaN OUTSIDE THE MASK
for k=1:num_images
    for i=1:(size(ellipse{k},3)-1);
        res{k,i}=rings{k}(:,:,i).*int_final{k};
        for l=1:size(int_final{k},1)
            for j=1:size(int_final{k},2)
            
                if rings{k}(l,j,i)==0;
                res{k,i}(l,j)=NaN;
                end
            end
        end
        mean{k,i}=nanmean(nanmean((res{k,i})));
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Test a possible gradient as a function of the distance from the spindle 
x_test=1:1:7;
for i=1:7
    y_test(i)=mean{1,i};
end
figure
plot(x_test,y_test);    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


