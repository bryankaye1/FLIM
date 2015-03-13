function [IllProfCal] = FLIM_IllProfCal(file)
% Input file of sdt acquisition of NADH or FAD. This smoothes the file to 
% get an illumination profile matrix that can be used to scale the eggs. 
% Important for calculating correct irradiance values for different eggs.
% file = 'Z:\Tim\2014-09-11,12 Cemb Mutants 2 no stripes\2014-09-15\IllProfCal\';

sdt = bh_readsetup([file]);
field = bh_getdatablock(sdt,1);
field = squeeze(sum(field,1));
% % for multiple files in a folder... obsolete
% if (isdir(file)&file(end)~='\') file = [file '\']; end
% D=dir([file '\*.sdt']); L = size(D,1);
% block_1=1;
% block_2=2;
% sdt = bh_readsetup([file D(1).name]);
% field = zeros(sdt.SP_IMG_Y,sdt.SP_IMG_X);
% for i = 1:L
%     sdt = bh_readsetup([file D(i).name]);
%     ch_1 = bh_getdatablock(sdt,block_1);
%     ch_1 = squeeze(sum(ch_1,1));
%     field = field + double(ch_1./L);
% end


%% Load non-uniform illumination field calibration file
n = round(size(field,1)/30);
%G = fspecial('gaussian',[n n],50);
close all;
%% smim = smooth2(field,50,50);
avefield2 = [field(:,n:-1:1) field field(:,end:-1:end-n+1)];
avefield2 = [avefield2(n:-1:1,:); avefield2; avefield2(end:-1:end-n+1,:)];
%% gim = imfilter(double(avefield2),G,'same');
gim = bpassTS(double(avefield2),n/2,512);
IllProfCal = gim(n+1:end-n,n+1:end-n);
IllProfCal = double(IllProfCal);
%% gim = bpassTS(field,10,512);
%figure;imshow(field,[]);
figure;imshow(IllProfCal,[]);%figure;imshow(smim,[]);
%save([UpOneDir(file) '\IllProfCal.mat'],'IllProfCal');

end