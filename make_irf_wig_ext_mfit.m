%[bneed, pulsewb, irf, irfsim, irf_pdf, wig, tmini, tmaxi, ext,varargout] =...
%    make_irf_wig_ext_minres(pth_irf, irfname, pth_wigs, wigsname, pth_ext, extname, data_shift_name, pth_data_for_shift,...
%    sfrc, shiftmin, shiftmax, w2step, w2min, w2max, backstep, backmin, backmax)

%This file is used to make the IRF vector, wiggles vector, extract vector. It also calculates the
%amount of missing time. It also makes the IRF_pdf (which is used for
%simulating data).


irfname = 'irf';
data_shift_name = 'atto565-10um';
wigsname = 'wigs';
cpath = '/Users/bryankaye/Documents/MATLAB/data/2015-5-21/';

pth_irf= cpath;
pth_wigs = cpath;
pth_data_for_shift = cpath;

sfrc = .2;

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

%% IRF and Wig shift
%Checks if files passed is in time series, if so, uses loads whole series
tsm = length(dir([pth_data_for_shift,data_shift_name,'_c*']));
[data1,~,~] = spc_2_his(tmini,tmaxi,data_shift_name,pth_data_for_shift,1,tsm);
data = data1/sum(data1);

intx = 1:sfrc:bins;
irfint = interp1(1:length(irf),irf,intx);
binskeep = bins - bneed;
s = Tlaser/bins:Tlaser/bins:Tlaser;
x = 1:binskeep;

p0 =[0,3.7,.02];
shifti= 0;
w2i = 3.7;
backi = 0.01;
[y] = flim_decay_shift(shifti, w2i,backi,sfrc, bins, bneed, irf, Tlaser, wig, data);



ft = fittype( 'flim_decay_shift(shifti, w2i, backi, sfrc, bins, bneed, irf, Tlaser, wig, x)',...
    'problem',{'sfrc', 'bins', 'bneed', 'irf', 'Tlaser', 'wig'},'independent','x' );

fit(x', data', ft, 'StartPoint', p0, 'problem',{sfrc, bins, bneed, irf, Tlaser, wig});

%Add a line that saves
%shiftb = shifti;
%w2b
%backb
%modelb
%irf
   
    %% Plots & irf 
    fprintf('shift %3.4f  lifetime %3.3f  back %3.6f Kres %3.6f \n', shiftb, w2b, backb, 1000*sum(res.^2));
    figure(11); clf; 
    subplot(2,1,1); hold on; plot(data); plot(modelb, 'r'); title('Shift: Model vs Data');
    subplot(2,1,2); plot(resb); title('Residue'); drawnow;


