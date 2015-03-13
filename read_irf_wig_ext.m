function [bneed, pulsewb, irf, irfsim, irf_pdf, wig, tmini, tmaxi, ext] = read_irf_wig_ext
%This file is used to make the IRF vector, wiggles vector, extract vector. It also calculates the
%amount of missing time. It also makes the IRF_pdf (which is used for
%simulating data). 

%Make sure the first 25 bins are a good representation of the noise for the IRF


Tlaser = 12.58;
%%%%%%% Read in IRF %%%%%%%%%
pth_sdt = 'Z:\bkaye\Bayes_2013\data\2014-1-22\'; %Type the file path of the Instrument Response Function (IRF) here.
file_name = 'irf-335pm.sdt'; irfname = file_name; %type the file name of the IRF here
data=dir([pth_sdt file_name]);
block=1; 
sdt = bh_readsetup([pth_sdt data.name]);
ch = bh_getdatablock(sdt,block); %raw lifetime data
ld = ch; %raw lifetime data in single form (normally stored as 4 bit number)
ld = double(ld);

%%%%%%%%% find the bins in which the system can record photons %%%%%%
bins = length(ld); %total number of bins
tmini = find(ld, 1 );
tmaxi = find(ld, 1, 'last' );
irf1 = ld(tmini:tmaxi);

%%%%%%%% calculate the number of bins you are missing %%%%%%%%
Tgraph = 14; %Recording interval of system (you may have to change the 14ns if you change the FLIM recording properties)
tpb = Tgraph/bins; %time width of one bin (time per bin)
pulsewb = round(Tlaser/tpb); %Number of bins corresponding to one period of the laser (12.58ns for one period of the laser)
bneed = pulsewb - length(irf1); %Number of bins we need to add to make one period

irf2 = irf1; 
    if bneed>0
        irf2 =  [irf2; zeros(bneed,1)]; %add bins (with 0 in each added bin) to make up 1 period
    end
  

%% Wig
%%%%%%%%% Read in Wigs %%%%%%%%%%%
pth_sdt = 'Z:\bkaye\Bayes_2013\data\2013-11-6\';
file_name = 'wigs.sdt';
data=dir([pth_sdt file_name]);
block=1; 
sdt = bh_readsetup([pth_sdt data.name]);
ch = bh_getdatablock(sdt,block); %raw lifetime data
ld = ch; %raw lifetime data in single form (normally stored as 4 bit number)
ld = double(ld);
wig1 = ld(tmini:tmaxi);
wig1 = wig1/mean(wig1); % This wigs get sent to bayes code to be put in the likelihood
m = 8;
wig =  tsmovavg(wig1','s',m)';
wig = [wig1(1:m-1); wig(m:end)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wig2 = wig; % Wig2 is the function that divides the wiggles out of the irf
    if bneed>0
        wig2 =  [wig2; ones(bneed,1)];  %add bins (with 1 in each added bin) to make up 1 period
    end


%%%%%%%%%%%%  Account for background noise and wigs in IRF  %%%%%%%%%%
irf2 = irf2-mean(irf2(1:25)); %calculates the background noise by taking the mean of the first 
irf2(irf2<0) =0; %Set all negative values to 0
irf = irf2./wig2; % Divide out wiggles
irfsim = floor(irf); %irfsim needs integer numbers to build pdf

    
%% irf_pdf
[irf_pdf] = pdf_builder(Tlaser, irfsim);

%% ext
%%%This section reads in the background extract signal, divides it out by
%%%the wiggles, then does a running average to smooth out poisson noise.

ne = 3;%number of extract measurements
for i=1:ne
pth_sdt = 'Z:\bkaye\Bayes_2013\data\2013-5-16\';
file_name = sprintf('m13_scan%d.sdt',i); %name of file, note that each extract measurement is called m13_scan1, m13_scan2, etc)
data=dir([pth_sdt file_name]);
block=1; 
sdt = bh_readsetup([pth_sdt data.name]);
ch = bh_getdatablock(sdt,block); %raw lifetime data
ld = ch; %raw lifetime data in single form (normally stored as 4 bit number)
ld = double(ld);
ext1 = ld(tmini:tmaxi);
ext1 =ext1./wig;% Account for wiggles
%figure(i); clf; plot(log10(ext1)); title(sprintf('Ext%d',i));
ext2(i,:) = ext1;
end

ext2 = sum(ext2);
figure(1); clf; plot(log10(ext2)); title('Ext Sum');
ext=ext2;
n = 8; %number of bins to be averaged in moving average
ext =  tsmovavg(ext,'s',n); %creates the moving average vector, which is smaller than the real vector. 
ext = [ext2(1:n-1), ext(n:end)]; %moving average vector will be smaller than real extract vector, so we add back the first few missing time spots
figure(2); clf; plot(log10(ext)); title('Ext Ave');
%figure(1); clf; plot((ext)); title('Non-log Ext Ave');


 %%%%%%  subtract background from extract pdf  %%%%%% 
ext = ext-mean(ext(1:25)); %Visually inspect ext on log scale, estimate background amount, and subtract this from ext
ext(ext<0) =0; %Set all negative values to 0
%ext = round(ext);
figure(3); clf; plot(log10(ext)); title('Ext thresh');
figure(4); clf; plot((ext)); title('non-log Ext thresh');
ext = ext/sum(ext); %Normalize extract pdf

save('Z:\bkaye\Bayes_2013\mat_files\sysin-2014-1-22-335pm.mat', 'bneed', 'pulsewb', 'irf', 'irfsim', 'irf_pdf', 'wig', 'tmini', 'tmaxi', 'ext','irfname');




%% Plots & irf

figure(1); clf; plot(log10(irf1)); title('IRF before threshold');
figure(2); clf; plot(log10(ext)); title('Ext');
figure(3); clf; plot(wig); title('wig');
figure(4); clf; plot(log10(irf)); title('IRF after threshold and de-wiggled');

end