function [bneed, pulsewb, irf, irfsim, irf_pdf, wig, tmini, tmaxi, ext] =...
    make_irf_wig_ext(pth_irf, irfname, pth_wigs, wigsname, pth_ext,...
    extname, data_shift_name, pth_data_for_shift)
%This file is used to make the IRF vector, wiggles vector, extract vector. It also calculates the
%amount of missing time. It also makes the IRF_pdf (which is used for
%simulating data).

%Make sure the first 25 bins are a good representation of the noise for the IRF

%irf_pdf removed
TOTALbins = 4096; %total number of bins
Tlaser = 12.58;
%%%%%%% Read in IRF %%%%%%%%%
[ld,~,~] =spc_2_his(1,4096,irfname,pth_irf,1,1);

%%%%%%%%% find the bins in which the system can record photons %%%%%%
tmini = find(ld, 1 );
tmaxi = find(ld, 1, 'last');
irf1 = ld(tmini:tmaxi)';
%%%%%%%% calculate the number of bins you are missing %%%%%%%%
Tgraph = 14; %Recording interval of system (you may have to change the 14ns if you change the FLIM recording properties)
tpb = Tgraph/TOTALbins; %time width of one bin (time per bin)
bins = round(Tlaser/tpb); %Number of bins corresponding to one period of the laser (12.58ns for one period of the laser)
bneed = bins - length(irf1); %Number of bins we need to add to make one period

irf2 = irf1;
if bneed>0
    irf2 =  [irf2; zeros(bneed,1)]; %add bins (with 0 in each added bin) to make up 1 period
end

%% Wig
%%%%%%%%% Read in Wigs %%%%%%%%%%%
[wig0,~,~] =spc_2_his(tmini,tmaxi, wigsname,pth_wigs,1,1);
wig1 = wig0'/mean(wig0); % This wigs get sent to bayes code to be put in the likelihood

if length(wig0) < 257
    wig = wig1;
else
    m = 8;
    wig =  tsmovavg(wig1','s',m)';
    wig = [wig1(1:m-1); wig(m:end)];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wig2 = wig; % Wig2 is the function that divides the wiggles out of the irf
if bneed>0
    wig2 =  [wig2; ones(bneed,1)];  %add bins (with 1 in each added bin) to make up 1 period
end
%%%%%%%%%%%%  Account for background noise and wigs in IRF  %%%%%%%%%%
irf2 = irf2-mean(irf2(1:round(length(irf2)/100))); %calculates the background noise by taking the mean of the first 1% of irf signal
irf2(irf2<0) =0; %Set all negative values to 0
%irf2(1500:end)=0; %This line removes irf past a certain time to remove false reflections etc
irf = irf2./wig2; % Divide out wiggles
irfsim = floor(irf); %irfsim needs integer numbers to build pdf


%% irf_pdf
%[irf_pdf] = pdf_builder(Tlaser, irfsim);
irf_pdf = 0;

%% ext
%%%This section reads in the background extract signal, divides it out by
%%%the wiggles, then does a running average to smooth out poisson noise.

%[ext1,~,~] =spc_2_his(tmini,tmaxi,extname,pth_ext, 1,1);
%ext2 =(ext1'./wig)';% Account for wiggles

%ext=ext2;
%n = 8; %number of bins to be averaged in moving average
%ext =  tsmovavg(ext,'s',n); %creates the moving average vector, which is smaller than the real vector.
%ext = [ext2(1:n-1), ext(n:end)]; %moving average vector will be smaller than real extract vector, so we add back the first few missing time spots

%%%%%%  subtract background from extract pdf  %%%%%%
%ext = ext-mean(ext(1:25)); %Visually inspect ext on log scale, estimate background amount, and subtract this from ext
%ext(ext<0) =0; %Set all negative values to 0
%ext = ext/sum(ext); %Normalize extract pdf
%ext1 = 1;


%% IRF and Wig shift
[data1,~,~] = spc_2_his(tmini,tmaxi,data_shift_name,pth_data_for_shift,1,1);
data1 = data1-mean(data1(1:10));
irf = irf - mean(irf(1:10));

data = data1/max(data1);
ibeg = find(data>.03,1,'first');
iend = find(data>.06,1,'first');

irf3 = irf/max(irf);
ibegirf = find(irf3>.03,1,'first');
iendirf = find(irf3>.06,1,'first');

shift1 = ibeg-ibegirf;
shift2 = iend-iendirf;
shiftb = round((shift1+0.2*shift2)/1.2);
irfb = circshift(irf3,[shiftb,0]);
gab = circshift(irf,[shiftb,0]);

%% Plots & irf

fprintf('shift %3.4f \n', shiftb);
%figure(13); clf; plot((data1),'b'); hold on; plot((gab), 'r'); title('Data (blue) & shifted IRF (RED)');
%figure(11); clf; plot(data(ibeg:iend),'b'); hold on; plot(irf3(ibeg:iend), 'r'); title('IRF and Data Sections');
figure(12); clf; plot((data),'b'); hold on; plot((irfb), 'r'); title('Normalized: Data (blue) & after shift IRF (RED)');
figure(13); clf; plot((data),'b'); hold on; plot(irf/max(irf), 'r'); title('Normalized: Data (blue) & before shift IRF (RED)');

pulsewb = bins;
wigsb = wig';
ext = wig';
save(strcat(pth_irf,'current.mat'), 'bneed', 'pulsewb', 'irf', 'irfsim', 'irf_pdf', 'tmini', 'tmaxi', 'ext','irfname','wigsb','gab');

end