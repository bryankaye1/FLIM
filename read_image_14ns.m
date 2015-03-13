

 pth_irf = 'C:\Users\Bryan\Documents\MATLAB\Bayes_2014\data\2014-9-5\';

data=dir([pth_sdt file_name]); %does this line do anything?
block=1; 
sdt = bh_readsetup([pth_sdt data.name]);
ch = bh_getdatablock(sdt,block); %raw lifetime data
ld = ch; %raw lifetime data in single form (normally stored as 4 bit number)
ld = double(ld);

sild = size(ld); %checks to see if this is an image (3dimensional) or just a FIFO (2 dimensional)
if length(sild) == 3
    ld = squeeze(sum(sum(ld,3),2)); %Averages all the pixels together
end

data = ld(tmini:tmaxi);
data=data';


% figure(1); clf; plot(x,ld); hold on; plot(tmin,ld(tmini),...
%     '--rs','MarkerSize',5); hold on; plot(tmax,ld(tmaxi), '--rs','MarkerSize',5);
% 
% figure(2); clf; plot(x,log10(ld));hold on; plot(tmin,log10(ld(tmini)),...
%     '--rs','MarkerSize',5); hold on; plot(tmax,log10(ld(tmaxi)), '--rs','MarkerSize',5);
% figure(3); clf; plot(irf2);
% figure(4); clf; plot(log10(irf2));

end