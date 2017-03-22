function [masks,mask_distance,angle_offset] = ellips_seg(int_image0,...
    pixel_length,scan_mag,rot_method)

mean_im0 = mean2(int_image0);
std_im0 = std2(int_image0);
thresh = mean_im0+(4/scan_mag)*std_im0;%to do: use a function to make the mask so 
%changes in registration and segmentation are gauraenteed to be the same
%IMPLEMENT make_masks if you want to use this function

thr_mask = int_image0>thresh;

%apply erosion and then dilation
mi0_pre = imfill(thr_mask,'holes');
se_erode = strel('disk',scan_mag/2,0);
er_mask = imerode(mi0_pre,se_erode);
se_dialate = strel('disk',scan_mag/2,0);
mi0 = imdilate(er_mask,se_dialate);
mi0_props = regionprops(mi0,'centroid','MajorAxisLength',...
    'MinorAxisLength','Orientation');
centroid(1) = mi0_props.Centroid(1)-65;
centroid(2) = mi0_props.Centroid(2)-65;

mi = mi0;
int_image = int_image0;

angle_offset=mi0_props.Orientation;
%%%%%%%%%%%%%%%%%%% in order to check this offset - create overlay first
%%%%%%%%%%%%%%%%%%% image, elliposid, MajorAxis and horizontal image axis;
phi_t = linspace(0,2*pi,50);
cosphi_t = cos(phi_t);
sinphi_t = sin(phi_t);


xbar_t = mi0_props.Centroid(1);
ybar_t = mi0_props.Centroid(2);

a_t = mi0_props.MajorAxisLength/2;
b_t = mi0_props.MinorAxisLength/2;

theta_t = pi*mi0_props.Orientation/180;
R_t = [ cos(theta_t)   sin(theta_t)
    -sin(theta_t)   cos(theta_t)];

xy_t = [a_t*cosphi_t; b_t*sinphi_t];
xy_t = R_t*xy_t;

x_t = xy_t(1,:) + xbar_t;
y_t = xy_t(2,:) + ybar_t;

hlen = mi0_props.MajorAxisLength/2;
cosOrient_tt = cosd(mi0_props.Orientation);
sinOrient_tt = sind(mi0_props.Orientation);
xcoords = xbar_t + hlen * [cosOrient_tt -cosOrient_tt];
ycoords = ybar_t + hlen * [-sinOrient_tt sinOrient_tt];
%line(xcoords, ycoords);
x_ax=round(mi0_props.Centroid-mi0_props.MajorAxisLength):1:...
    round(mi0_props.Centroid+mi0_props.MajorAxisLength);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%% VALUES to rotate: theta_max + angle_offset from first image!
int_final=imrotate(int_image ,(-angle_offset),rot_method,'crop');
mi_final=imrotate(mi,(-angle_offset),rot_method,'crop');
%%%%% HERE: all spindles are centered and aligned along the image axis;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% next: create ring masks
%mag_frac=.7:0.05:3;
i = 0;
kvec = -round(.3*mi0_props.MinorAxisLength/2):64;
for k=kvec
    i = i+1;
    %We want to expand the elipse in units of single pixels along the short
    %axis. If the short, long axis is a,b, then expansions are set such
    %that a/b = (a+k) / (b+ n ), where k is an integer. This forces n to be
    %k*b/a
    a = mi0_props.MinorAxisLength/2;
    b = mi0_props.MajorAxisLength/2;
    new_a(i) = mi0_props.MinorAxisLength/2 + k;
    new_b(i) = mi0_props.MinorAxisLength/2 + k;
    
    %%% only allow for ellipses that fit the image:
    if new_b(i) > 64
    else
        %%% create masks of the respective ellipses for each images:
        [X Y]= meshgrid(1:128,1:128);
        % ellipse(:,:,i) = ((X-65)/(elip_len(i)*mi0_props.MajorAxisLength/2)).^2 + ...
        %     ((Y-65)/(elip_len(i)*mi0_props.MinorAxisLength/2)).^2<=1;
        ellipse(:,:,i) = ( (X-65)/new_b(i) ).^2 + ( (Y-65)/new_a(i) ).^2 <= 1;
    end
end

%%% create masks:
for i=1:(size(ellipse,3)-1)
    masks (:,:,i) = ellipse(:,:,i+1)-ellipse(:,:,i);
    if kvec(i) == 0
        spindle_edge = find(masks(65,:,i),1);
    end
end

for i=1:(size(ellipse,3)-1)
    mask_distance(i) = pixel_length*(spindle_edge - find(masks(65,:,i),1));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% ATTEMPT TO VISUALIZE POSITION OF masks ON TOP OF THE INTENSITY

h3 = figure(3); clf;
B= imoverlay(mat2gray(int_final),max(masks(:,:,1),masks(:,:,end)));
imshow(B,'InitialMagnification','fit');
title('aligned registered image and first and last mask');
end
% for i=[1,(size(ellipse,3)-1)]
%     
%     R_norm = [ 1,0;0,1];
%     
%     xy_norm = [new_b(i)*cosphi_t; new_a(i)*sinphi_t];
%     xy_norm = R_norm*xy_norm;
%     
%     x_norm(:,:,i) = xy_norm(1,:) + 65;
%     y_norm(:,:,i) = xy_norm(2,:) + 65;
%     
%     plot(x_norm(:,:,i),y_norm(:,:,i),'r','LineWidth',2);
%     title('aligned registered image and first and last mask');
% end



%%%%%%%%%%%%%%%%%%%%%%% OVERLAY THESE RINGS WITH IMAGES;
%%%%%%%%%%% IN CASE OF MEAN - SET "res" TO NaN OUTSIDE THE MASK
% for i=1:(size(ellipse,3)-1)
%     res{i}=rings(:,:,i).*int_final;
%     for l=1:size(int_final,1)
%         for j=1:size(int_final,2)
%
%             if rings(l,j,i)==0
%                 res{i}(l,j)=NaN;
%             end
%         end
%     end
%     mean_res{i}=nanmean(nanmean((res{i})));
% end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%% Test a possible gradient as a function of the distance from the spindle
% x_test=1:1:7;
% for i=1:7
%     y_test(i)=mean{1,i};
% end
% figure
% plot(x_test,y_test);