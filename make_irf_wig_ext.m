function [bneed, pulsewb, irf, irfsim, irf_pdf, wig, tmini, tmaxi, ext] =...
    make_irf_wig_ext(pth_irf, irfname, pth_wigs, wigsname, pth_ext, extname, data_shift_name, pth_data_for_shift,...
    sfrac, shiftmin, shiftmax, w2step, w2min, w2max, backstep, backmin, backmax, wigstep, wigmin, wigmax)
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
[ext1,~,~] =spc_2_his(tmini,tmaxi,extname,pth_ext, 1,1);
ext2 =(ext1'./wig)';% Account for wiggles

ext=ext2;
n = 8; %number of bins to be averaged in moving average
ext =  tsmovavg(ext,'s',n); %creates the moving average vector, which is smaller than the real vector.
ext = [ext2(1:n-1), ext(n:end)]; %moving average vector will be smaller than real extract vector, so we add back the first few missing time spots

%%%%%%  subtract background from extract pdf  %%%%%%
ext = ext-mean(ext(1:25)); %Visually inspect ext on log scale, estimate background amount, and subtract this from ext
ext(ext<0) =0; %Set all negative values to 0
ext = ext/sum(ext); %Normalize extract pdf

%% IRF and Wig shift

if strcmp(data_shift_name,'manual')  

   tempf =  load('ll-his');
   dtemp = tempf.info.his;
   data1 = dtemp(tmini:tmaxi)';
   
%     ndiv = 1;
%     [data1, ~] = read_spc_timeseries([pth_data_for_shift data_shift_name],ndiv);
%     [data1] =   data1(tmini:tmaxi);
else
    [data1,~,~] = spc_2_his(tmini,tmaxi,data_shift_name,pth_data_for_shift,1,1);
end

intx = 1:sfrac:bins;
irfint = interp1(1:length(irf),irf,intx);
sumres = inf;
binskeep = bins - bneed;

for wigshifti = wigmin:wigstep:wigmax
    
    wigs = circshift(wig,wigshifti)';
    
    for shifti = shiftmin:sfrac:shiftmax
        
        shift = round(shifti/sfrac);
        irf2 = circshift(irfint',shift);
        ga = interp1(intx, irf2, 1:bins);
        
        for w2i = w2min:w2step:w2max
            
            s = Tlaser/bins:Tlaser/bins:Tlaser;
            f2 = exp(-s/w2i); %signal over one period
            f2 = [f2 f2]; %signal over 2 consecutive periods
            f2con = conv(f2,ga); %PDF after conv
            f2bar = f2con(bins+1:2*bins); %pdf after mod-ing by laser period
            f2h = f2bar(1:binskeep); %Keep only the appropriate bins
            f2h = f2h/sum(f2h);
            for backi = backmin:backstep:backmax
                back = backi/bins;
                model = (f2h + back).*wigs;
                
                model = model/sum(model);
                data = data1/sum(data1);
                ga =ga/sum(ga);
                
                res = (data-model)./sqrt(data);
                sumresi = sum(res.^2);
                
                if sumresi < sumres
                    sumres = sumresi;
                    shiftb = shifti;
                    gab = ga;
                    w2b = w2i;
                    backb = backi;
                    wigshiftb = wigshifti;
                    wigsb = wigs;
                    
                    modelb = model;
                    datab = data;
                    resb = res;
                    
                end
            end
        end
    end 
    
end
    
    %% Plots & irf
    
    fprintf('shift %3.4f  lifetime %3.3f  back %3.6f wigshift %3.1f Kres %3.6f \n', shiftb, w2b, backb, wigshiftb, 1000*sum(res.^2));
    figure(11); clf; plot((datab)); hold on; plot((modelb), 'r'); title('Shift: Model vs Data');
    figure(12); clf; plot(resb); title('Residue');
    
%     figure(3); clf; plot(log10(irf1)); title('IRF before threshold');
%     figure(4); clf; plot(log10(irf)); title('IRF after threshold and de-wiggled');
%     figure(5); clf; plot(wig); title('wig');
%     figure(6); clf; plot(log10(ext)); title('Ext');
%     
    pulsewb = bins;
    save('Z:\bkaye\Bayes_2013\mat_files\current.mat', 'bneed', 'pulsewb', 'irf', 'irfsim', 'irf_pdf', 'tmini', 'tmaxi', 'ext','irfname','wigsb','gab');
    
end