
close all;
pth_data1 = 'C:\Users\Bryan\Documents\MATLAB\data\2015-2-14\';
pth_data2 = 'C:\Users\Bryan\Documents\MATLAB\data\2015-2-25\';
dataname1 = 'e1-m1-s1';
dataname2 = 'e1-m5-s1';

[pmat,jmax,ni,sinti,imagedata1] = spc_2_his(1,4096,dataname1,pth_data1,1,1,'imout');
im1 = squeeze(imagedata1);

[pmat,jmax,ni,sinti,imagedata2] = spc_2_his(1,4096,dataname2,pth_data2,1,1,'imout');
im2 = squeeze(imagedata2);
%%
whiteval = max([max(max(im1)) max(max(im2))]);
blackval = min([min(min(im1)) min(min(im2))]);

I1 = mat2gray(im1,[blackval whiteval]);
I2 = mat2gray(im2,[blackval whiteval]);

figure(1); clf; imshow(I1); title(dataname1);
figure(2); clf; imshow(I2); title(dataname2);


%For Pf figure, dataname1 = 'm2_c20' - dataname2 = 'm4_c20';
% For nPol figure, 