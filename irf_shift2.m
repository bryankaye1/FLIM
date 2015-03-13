%% Load-in

tempf=load('Z:\bkaye\Bayes_2013\mat_files\sysin-2014-1-22-335pm.mat','-mat','irf', 'bneed', 'pulsewb','wig','tmini','tmaxi','ext');
irf=tempf(1).irf;        %Probability of delay (pdf_IRF)
brem = tempf(1).bneed;   %Number of bins to remove from the 12.5ns to match BH data
bins = tempf(1).pulsewb; %Number of bins that make up 12.5ns
wig = (tempf(1).wig)';   %wig = ones(1,length(wig)); fprintf('wigs set to 0 in fit but not data\n');
tmini = tempf(1).tmini;
tmaxi = tempf(1).tmaxi;
ext = tempf(1).ext;
binskeep = bins-brem;

% pth_sdt = 'Z:\bkaye\Bayes_2013\data\2013-10-11/';
% file_name = 'coum_2e5_p1pc_m1.asc';%'erthy_2e5_p25pc_m1.asc';
% [hislong,~, ~] = read_spc_fn([pth_sdt file_name],0);

%data = hislong(tmini:tmaxi)';

pth_sdt = 'Z:\bkaye\Bayes_2013\data\2014-1-22/';
file_name = 'm0.sdt';
[data] = read_data_14ns(tmini,tmaxi,file_name,pth_sdt);%p is the histogram data. file_name is the name of the data file



%% Main Code


sfrac = .1;
intx = 1:sfrac:bins;
irfint = interp1(1:bins,irf,intx);
sumres = inf;
binskeep = pulsewb - bneed;


for shifti = -8:sfrac:-5
    
    shift = round(shifti/sfrac);
    irf2 = circshift(irfint',shift);
    ga = interp1(intx, irf2, 1:bins);
    
    wigshift = -15;
    wigs = circshift(wig',wigshift)';
    T = 12.58;
    
    for w2i = 2.5:.02:4
        %for backi = 0:.0001:.001
        %for wigshift = -10:1:-4
        
        s = T/bins:T/bins:T;
        f2 = exp(-s/w2i); %signal over one period
        f2 = [f2 f2]; %signal over 2 consecutive periods
        f2con = conv(f2,ga); %PDF after conv
        f2bar = f2con(bins+1:2*bins); %pdf after mod-ing by laser period
        f2h = f2bar(1:binskeep); %Keep only the appropriate bins
        f2h = f2h/sum(f2h);
        for backi = .05:0.005:.1
            back = backi/bins;
            model = (f2h + back).*wigs;
            
            model = model/sum(model);
            data = data/sum(data);
            ga =ga/sum(ga);
            
            res = (data-model)./sqrt(data);
            sumresi = sum(res.^2);
            
            if sumresi < sumres
                sumres = sumresi;
                shiftb = shifti;
                w2b = w2i;
                backb = backi;
                wigshiftb = wigshift;
                
                modelb = model;
                datab = data;
                resb = res;
                
            end
        end
        %end
    end
end
fprintf('shift %3.4f  lifetime %3.3f  back %3.6f wigshift %3.1f Sres %3.6f \n', shiftb, w2b, backb, wigshiftb, sum(abs(resb)));
figure(1); clf; plot((datab)); hold on; plot((modelb), 'r');% hold on; plot(log(ga),'g');
figure(2); clf; plot(resb);