function [data] = read_data_14ns(tmini,tmaxi,file_name,pth_sdt)

data=dir([pth_sdt file_name]); %does this line do anything?
block=1; 
sdt = bh_readsetup([pth_sdt data.name]);
ch = bh_getdatablock(sdt,block); %raw lifetime data
ld = ch; %raw lifetime data in single form (normally stored as 4 bit number)
ld = double(ld);

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