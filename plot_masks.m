function plot_masks(seg_results,varargin)
%add functions instructions here 

image_stack = seg_results.offset_corrected_image_stack;
masks = seg_results.int_masks;

if nargin == 2
    close_mask = varargin{1};
elseif nargin == 3
    close_mask = varargin{1};
    far_mask = varargin{2};
else
    close_mask = 1;
    [a,b,far_mask] = size(masks);
end

figure; clf;

B= imoverlay(mat2gray(image_stack),max( masks(:,:,close_mask),masks(:,:,far_mask )));
imshow(B,'InitialMagnification','fit');
title('"Equi-distance" masks');

end