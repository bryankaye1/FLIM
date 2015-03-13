%function [hislong,hisshort, data,dshortt] = read_spc_fn_v2(filename,mexp)
clear
%% Get Short-Lifetime Photon Master List
tic
%filename = 'Z:\bkaye\Bayes_2014\data\2014-10-7\coumarin_1p4e5_550-88nm_950ex_2000sec_m1.asc';
fileID = fopen('C:\Users\Bryan\Documents\MATLAB\Bayes_2014\data\2014-10-7\eryth-4p6e6_550-88nm-950ex-100sec_m1.asc');
data1 = textscan(fileID, '%f %f %f %f','HeaderLines', 19 );
fclose(fileID);

dshort = data1{2}; %Long lifetime photons (arrival times): dlong = distribution long
npshort = length(dshort); %number of photons
toc

edges = 0:4095;
hisshort = histc(dshort, edges);

save('sh-his','hisshort');
save('sh-at','dshort');

%% Get Long-Lifetime Photon Master List%%
tic
fileID = fopen('C:\Users\Bryan\Documents\MATLAB\Bayes_2014\data\2014-10-7\coumarin_4p3e6_550-88nm_950ex_200sec_m1.asc');
data2 = textscan(fileID, '%f %f %f %f','HeaderLines', 19 );
fclose(fileID);
dlong = data2{2};
nplong = length(dlong);
toc

edges = 0:4095;
hislong = histc(dlong, edges);

save('lh-his','hislong');

save('lh-at','dlong');