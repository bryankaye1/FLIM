function [bneed pulsewb irf irfsim irf_pdf wig tmini tmaxi] = read_irf_wig
%This file is used to read in the IRF and wiggles. It also calculates the
%amount of missing time.

%IRF for sim uses
%%%%%%% Read in IRF %%%%%%%%%
pth_sdt = 'C:\Users\Bryan\Documents\MATLAB\Bayes_2013\data\2013-4-18\';
file_name = 'irf1.sdt';
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
Tgraph = 14; %Recording interval of system
tpb = Tgraph/bins; %time width of one bin (time per bin)
pulsewb = round(12.5/tpb); %Number of bins corresponding to one period of the laser
bneed = pulsewb - length(irf1); %Number of bins we need to add to make one period

irf2 = irf1; 
    if bneed>0
        irf2 =  [irf2; zeros(bneed,1)]; %add bins (with 0 in each added bin) to make up 1 period
    end
  

%% Wig
%%%%%%%%% Read in Wigs %%%%%%%%%%%
pth_sdt = 'C:\Users\Bryan\Documents\MATLAB\Bayes_2013\data\2013-5-10\';
file_name = 'wiggles_300m.sdt';
data=dir([pth_sdt file_name]);
block=1; 
sdt = bh_readsetup([pth_sdt data.name]);
ch = bh_getdatablock(sdt,block); %raw lifetime data
ld = ch; %raw lifetime data in single form (normally stored as 4 bit number)
ld = double(ld);
wig1 = ld(tmini:tmaxi);
wig1 = wig1/mean(wig1); % This wigs get sent to bayes code to be pu tin the likelihood
m = 4;
wig =  tsmovavg(wig1','s',m)';
wig = [wig1(1:m-1); wig(m:end)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wig2 = wig; % Wig2 is the function that divides the wiggles out of the irf
    if bneed>0
        wig2 =  [wig2; ones(bneed,1)];  %add bins (with 1 in each added bin) to make up 1 period
    end


    %%%%  thresholding  %%%% 
irf2 = irf2-round(10^2.3); %Visually inspect IRF on log scale, estimate background amount, and subtract this from IRF
irf2(irf2<0) =0; %Set all negative values to 0
irf = irf2./wig2; % Divide out wiggles
irfsim = floor(irf); %irfsim needs integer numbers to build pdf

    
%% irf_pdf
[irf_pdf] = pdf_builder(12.5, irfsim);

%% Plots & irf

figure(1); clf; plot(log10(irf1)); title('IRF before threshold');
figure(2); clf; plot(log10(irf2)); title('IRF after threshold');
figure(5); clf; plot(wig); title('wig');
figure(4); clf; plot(log10(irf)); title('IRF after threshold and de-wiggled');

end