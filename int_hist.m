%% Clears and Comments
%Here we prepare the data and set the analysis parameters
%pr means FRET fraction now
%clear;
clear;


for el = 1:1
    jmax = 1;
    exptmax = 1;
    cyclesmax = 1;
    tfw = .6;%(el-1)*.1;
    tbac = .2;
    
    %jind =0;
    for jind = 1:jmax %while jind < jmax%
        % jind = jind +1;
        
        for expt = 1:exptmax
            if expt == 1
                   str(1).name = strcat('m3n-spot1-ll.sdt');
                
            end
            for cindex = 1:cyclesmax
                %% Load in IRF, wigs
                tempf=load('Z:\bkaye\Bayes_2013\mat_files\current.mat','-mat','bneed',...
                    'pulsewb', 'irf', 'irfsim', 'irf_pdf', 'tmini', 'tmaxi', 'ext','irfname','wigsb','gab');
                brem = tempf(1).bneed;   %Number of bins to remove from the 12.5ns to match BH data
                bins = tempf(1).pulsewb; %Number of bins that make up 12.5ns
                tmini = tempf(1).tmini;
                tmaxi = tempf(1).tmaxi;
                ext = tempf(1).ext;
                wig = tempf(1).wigsb;
                ga = tempf(1).gab;
                binskeep = bins-brem;
                
                pth_data = 'C:\Users\Bryan\Documents\MATLAB\Bayes_2014\data\2014-9-5\';
                dataname = str(cindex).name;
                
                data=dir([pth_data dataname]); %does this line do anything?
                block=1;
                sdt = bh_readsetup([pth_data data.name]);
                ch = bh_getdatablock(sdt,block); %raw lifetime data
                ld = ch; %raw lifetime data in single form (normally stored as 4 bit number)
                ld = double(ld);
                datatemp = ld(tmini:tmaxi,:,:);
                datat = ld;
                
                for j=1:128
                    for k = 1:128
                        pmt((j-1)*128+k) = sum(datat(:,j,k),1);
                    end
                end
               %% 
                
               xmin = min(pmt);
               xmax = max(pmt);
               xstep = (xmax - xmin)/100;
               
               xvec = round(xmin:xstep:xmax);
                
               [bincts, bind] = histc(pmt,xvec);
              % bincts = [0 bincts(2:end)];
                [ctsmax, mindex] = max(bincts);
                bmon = (mindex-1)*xstep+xmin;
                %bmon = 247;
                
                monc = poisspdf(xvec,bmon)*max(bincts)/max(poisspdf(xvec,bmon));
                %gfit = fit(xvec',bincts','gauss1');
              % monc = gfit.a1.*exp(-((xvec-gfit.b1)./gfit.c1).^2);
               polc = bincts - monc;              
                
              % figure(el); clf; plot(xvec,bincts); title(str(1).name);
%                 figure(11); clf; plot(xvec,bincts); title('total cts vs int');
%                figure(12); clf; plot(monc);  title('mon cts');
%                figure(13); clf; plot(xvec,monc);  title('mon cts vs int');
%                figure(14); clf; plot(polc); title('polc');
                
                 figure(3+el); clf; hold all; plot(bincts); plot(monc); plot(polc);title(str(1).name); legend('bincts','monc','polc');
                
%                  gev = [125,92,100,75,67,64]; 
%                  ge = gev(el); %ge(6) = 64; 67, 75,100,92,128
                 pho_pol = sum(polc.*xvec);
                 pho_mon = sum(monc.*xvec);
                 a = .3;
                 b_pol = 1-a*(1-1.5/3.87);
                 pol_frac = pho_pol/(pho_mon*b_pol+pho_pol);
                 fprintf('bmon is %1.2f\n',bmon);
                     pause(1);
                %figure(5); clf; plot(xvec,poisspdf(xvec,bmon));
               % figure(6);clf; plot(poisspdf(1:100,11));
                
            end
        end
    end
    
end

