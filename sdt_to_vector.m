function [ld] = sdt_to_vector(pth_sdt,file_name)

sdt = bh_readsetup([pth_sdt file_name]); block=1;
ch = bh_getdatablock(sdt,block); %raw lifetime data
ld = ch; %raw lifetime data in single form (normally stored as 4 bit number)
ld = double(ld); 
end