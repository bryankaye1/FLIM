function [mi0] = make_masks(int_image0,di_erode_pixels,thr_mult,back_vec,varargin)

if nargin ==4
image4thresh = int_image0;
image4thresh(image4thresh==0) = NaN;
thresh = nanmean(nanmean(image4thresh(back_vec,back_vec)))*thr_mult;
%thresh = mean2(int_image0(back_vec,back_vec)))*thr_mult;
else
image4thresh = varargin{1}{1};
image4thresh(image4thresh==0) = NaN;   
thresh = nanmean(nanmean(image4thresh(back_vec,back_vec)))*thr_mult;     
%thresh = mean2(varargin{1}{1}(back_vec,back_vec))*thr_mult; 
end
thr_mask = int_image0>thresh;
%apply erosion and then dilation
mi0_pre = imfill(thr_mask,'holes');
se_erode = strel('disk',di_erode_pixels,0);
er_mask = imerode(mi0_pre,se_erode);
se_dialate = strel('disk',di_erode_pixels,0);
mi0 = imdilate(er_mask,se_dialate);
mi0 = imfill(mi0,'holes');
end