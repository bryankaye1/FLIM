%% Clears and Comments
clear;
%clc;
clf;

%Inputs: Donor lifetime.
%Outputs: Background, FRET lifetime, FRET and Non-FRET fractions
%V8: added w1min
%Added extract fitting dimension.

%Need to add:
%1) For window factor, don't use best in each dimension, but instead use best in parameter space

%2) Marginalize to find frac_est, then apply scale factor
tic
for expt = 1:1
    
    if expt == 1
        str(1).name = 'm2_c1.sdt';
        str(2).name = 'm2_c2.sdt';
        str(3).name = 'm2_c01.sdt';
        
        str(4).name = 'm2_c02.sdt';
        str(5).name = 'm2_c03.sdt';
        str(6).name = 'm2_c04.sdt';
        
        str(7).name = 'm2_c05.sdt';
        str(8).name = 'm2_c06.sdt';
        str(9).name = 'm2_c07.sdt';
        
        str(10).name = 'm2_c08.sdt';
        str(11).name = 'm2_c09.sdt';
        str(12).name = 'm2_c10.sdt';
        
        str(13).name = 'm2_c11.sdt';
        str(14).name = 'm2_c12.sdt';
        str(15).name = 'm2_c13.sdt';
        
        str(16).name = 'm2_c14.sdt';
        str(17).name = 'm2_c15.sdt';
        str(18).name = 'm2_c16.sdt';
        
        str(19).name = 'm2_c17.sdt';
        str(20).name = 'm2_c18.sdt';
        str(21).name = 'm2_c19.sdt';
        
        str(22).name = 'm2_c20.sdt';
        str(23).name = 'm2_c21.sdt';
        str(24).name = 'm2_c22.sdt';

        str(25).name = 'm2_c23.sdt';
        str(26).name = 'm2_c24.sdt';
        str(27).name = 'm2_c25.sdt';
        
        str(28).name = 'm2_c26.sdt';
        str(29).name = 'm2_c27.sdt';
        str(30).name = 'm2_c28.sdt';        
        
        irfname = 'irf-1053pm.sdt';
        data_shift_name = 'm2_c1.sdt';
       
        w1min = .7; w1max = .7;
       w2min =  3.95; w2max =  3.95; % w2min =  3.8755; w2max =  3.8755;
        w3 = 0.6677;
        f3frac = 0.0199;
        
    else
        str(1).name = 'm0-e2-s1-spot1.sdt';
        str(2).name = 'm0-e2-s1-spot2.sdt';
        str(3).name = 'm0-e2-s1-spot3.sdt';
        
        str(4).name = 'm0-e2-s2-spot1.sdt';
        str(5).name = 'm0-e2-s2-spot2.sdt';
        str(6).name = 'm0-e2-s2-spot3.sdt';
        
        str(7).name = 'm0-e2-s3-spot1.sdt';
        str(8).name = 'm0-e2-s3-spot2.sdt';
        str(9).name = 'm0-e2-s3-spot3.sdt';
        
        irfname = 'irf-1121pm.sdt';
        data_shift_name = 'm0-e2-s2-spot2.sdt';
        
        w1min = .9487; w1max = .9487;
        w2min =  3.9050; w2max =  3.9050;  
        w3 =  0.6998;
        f3frac =  0.0201;
        
    end
    
    for index = 1:30






%% Make IRF, wigs, IRF shift, wigs shift, and extract signal
sysinfo = 0; % Set to 1 if you want to force a rerun of make_irf_wig_ext
pth_irf = 'C:\Users\Bryan\Documents\MATLAB\Bayes_2014\data\2014-4-20\'; %Type the file path of the Instrument Response Function (IRF) here.
%irfname = 'irf-641pm.sdt';

pth_wigs = 'Z:\bkaye\Bayes_2013\data\2014-2-13\';
wigsname = 'wigs_1p5e6.sdt'; 

pth_ext = 'Z:\bkaye\Bayes_2013\data\2013-5-16\';
extname = 'm13_scan1.sdt';

pth_data_for_shift = 'C:\Users\Bryan\Documents\MATLAB\Bayes_2014\data\2014-4-20\';
%data_shift_name = 'm0-e1-s1-spot1.sdt';

cppout = fopen('SysInfo.txt');
[old_filenames, count] = fscanf(cppout, '%s');
new_filenames = [pth_irf, irfname, pth_wigs, wigsname, pth_ext,extname, pth_data_for_shift, data_shift_name];

if isequal(new_filenames,old_filenames) && sysinfo == 0;
    
else
    sfrac = .2; shiftmin = -10; shiftmax = 0;
    w2step_shift = 0.01; w2min_shift = 3.5; w2max_shift = 4;
    backstep = 0.01; backmin = 0.03; backmax = 0.1;
    wigstep = 2; wigmin = 0; wigmax = 0;
    
    [~, ~, ~, ~, ~, ~, ~, ~, ~] = make_irf_wig_ext(pth_irf, irfname, pth_wigs, wigsname, pth_ext, extname, data_shift_name, pth_data_for_shift, sfrac, shiftmin, shiftmax, w2step_shift, w2min_shift, w2max_shift, backstep, backmin, backmax, wigstep, wigmin, wigmax);
    fileID = fopen('SysInfo.txt','w');
    fprintf(fileID,'%s%s\n%s%s\n%s%s\n%s%s\n',pth_irf, irfname, pth_wigs, wigsname, pth_ext,extname, pth_data_for_shift, data_shift_name);
    fclose(fileID);
    
end

%% Load in IRF, wigs

tempf=load('Z:\bkaye\Bayes_2013\mat_files\current.mat','-mat','irf','bneed',...
    'pulsewb', 'irf', 'irfsim', 'irf_pdf', 'tmini', 'tmaxi', 'ext','irfname','wigsb','gab');
irf=tempf(1).irf;        %Probability of delay (pdf_IRF)
brem = tempf(1).bneed;   %Number of bins to remove from the 12.5ns to match BH data
bins = tempf(1).pulsewb; %Number of bins that make up 12.5ns
tmini = tempf(1).tmini;
tmaxi = tempf(1).tmaxi;
ext = tempf(1).ext;
wig = tempf(1).wigsb;
ga = tempf(1).gab;
binskeep = bins-brem;

pth_data = 'C:\Users\Bryan\Documents\MATLAB\Bayes_2014\data\2014-4-20\';
dataname = str(index).name;
%dataname = 'm5-e1-s1-spot1.sdt'; ndiv = 1;

thrt = 0; %0 means do not threshold, 1 means threshold.
imthr = 3;%amount to threshol by. Thresh is max(int)/imthr
[p] = read_imagethr_14ns(tmini,tmaxi,dataname,pth_data,imthr,thrt);

% 
% [data, timeinterval] = read_spc_timeseries([pth_data dataname],ndiv);
% [p] =   data(tmini:tmaxi);  %hislong(tmini:tmaxi)';%
% 
% pth_data = 'Z:\bkaye\Bayes_2013\data\2014-2-3/';
% dataname = 'm0.sdt';
% [p] = read_data_14ns(tmini,tmaxi,dataname,pth_data);% p is the histogram data. dataname is the name of the data file


%%
for l1 = 1:2
    %% Scan Parameters and Simulation Variables/Load Data
    if l1 == 1
        w1step = .01; %w1min= 0.6; w1max =  1.3; %w1min must be an integer multiple of w1step.
        w2step = .01; %w2min =  3.8755; w2max =  3.8755;
        fracstep = 0.005;
        prstep = fracstep; prmin=0; prmax = 1; 
        w02step = fracstep; w02min = 0; w02max = 1;
        extstep = fracstep; extmin = 0; extmax = 0;
        thr = 0.01; 
        T = 12.58; %Used in both and sim and real data
        prstep0 = prstep;  w02step0 = w02step; w2step0 = w2step; w1step0 = w1step;
    else
        [w1step, w2step, prstep, w02step, extstep] =...
            newstepsize(w1max,w1min,w2max,w2min,prmax,prmin,w02max,w02min,extmax,extmin,lw1,lw2,lpr,lw02,lext);
    end
    
    %% Calculate Post
    
    %Initialize loglike
    loglike = -pi*ones(round(1+(w1max-w1min)/w1step),round(1+(w2max-w2min)/w2step),...
        round(1+(w02max-w02min)/w02step),...
        round(1+(prmax-prmin)/prstep),round(1+(extmax-extmin)/extstep));
    
    %Estimate time to search space
    sizell = size(loglike,1)*size(loglike,2)*size(loglike,3)*size(loglike,4)*size(loglike,5)/1000;
    time_est= ceil(.06*sizell);
    %fprintf('\n\n%1.0fk gridpoints. Runtime %1.0f seconds\n',sizell,time_est);
    
    s = T/bins:T/bins:T; %time vector used to generate PDF of signal exponential
    %GENERATE PDF FOR W2 (NON-FRET)
    %tic
    for w2i = w2min:w2step:w2max
        f2 = exp(-s/w2i); %signal over one period
        f2 = [f2 f2]; %signal over 2 consecutive periods
        f2con = conv(f2,ga); %PDF after conv
        f2bar = f2con(bins+1:2*bins); %pdf after mod-ing by laser period
        f2h = f2bar(1:binskeep); %Keep only the appropriate bins
        b = sum(f2h)/sum(f2bar); %Norm constant for entire period / Norm constant for window period
        f2h = f2h/sum(f2h); %Normalized PDF of detecting a photon in a bin given photon came from w2 (non-FRET)
             
        %w3 = 0.6677; % Found from e1 experiment. 
        f3 = exp(-s/w3); %signal over one period
        f3 = [f3 f3]; %signal over 2 consecutive periods
        f3con = conv(f3,ga); %PDF after conv
        f3bar = f3con(bins+1:2*bins); %pdf after mod-ing by laser period
        f3h = f3bar(1:binskeep); %Keep only the appropriate bins
        f3h = f3h/sum(f3h); %Normalized PDF of detecting a photon in a bin given photon came from w2 (non-FRET)    
        
        
        f2h = (f3frac* 1.0789)*f3h + f2h;
        f2h = f2h/sum(f2h);
        
        %SEARCH OVER W1,W02,W00
        for w1i = w1min:w1step:w1max
            %w1i
            clear f;
            f1 = exp(-s/w1i);
            f = [f1 f1];
            fcon = conv(f,ga);
            fbar = fcon(bins+1:2*bins);
            fh = fbar(1:binskeep);
            a = sum(fh)/sum(fbar); %Norm constant for entire period / Norm constant for window period
            fh = fh/sum(fh);
            %b = 1; a =1;
            prf = (a/b);%*(w1i/w2i);
            %search over all w02 and w00 space, but set impossible states
            %(w00 +w01 > 0) to zero probability (-inf in loglike)
            
            for w02i = w02min:w02step:w02max
                for pri = prmin:prstep:prmax
                    for exti = extmin:extstep:extmax
                        
                        if pri*prf*w02i+ w02i +exti > 1
                            loglike(round(1+(w1i-w1min)/w1step),round(1+(w2i-w2min)/w2step), round(1+(w02i-w02min)/w02step),...
                                round(1+(pri-prmin)/prstep),round(1+(exti-extmin)/extstep)) = -inf;
                        else
                            loglike(round(1+(w1i-w1min)/w1step),round(1+(w2i-w2min)/w2step), round(1+(w02i-w02min)/w02step),...
                                round(1+(pri-prmin)/prstep),round(1+(exti-extmin)/extstep)) =...
                                sum(log(wig.*((1-pri*prf*w02i-w02i-exti)/binskeep+(pri*prf*w02i)*fh + w02i*f2h + exti*ext)).*p);
                        end
                    end
                end
            end
        end
    end
    %toc
    
    checkloglike(loglike); %check's loglike for matrix size errors (makes sure entire space is scanned)
    loglike2 = loglike - max(max(max(max(max(loglike))))); %Sets max of like to 1.
    like = exp(loglike2);
    
    %%%%%% Priors %%%%%%
    pri = ones(round(1+(w1max-w1min)/w1step),round(1+(w2max-w2min)/w2step),round(1+(w02max-w02min)/w02step),round(1+(prmax-prmin)/prstep),round(1+(extmax-extmin)/extstep)); %Ini Properly
    post = like; %post = like.*pri; 
    
    
    %% Marg and Plots
    
    [prest, prBesti, prBest, prestx, w1est, w1Besti,w1Best, w1estx, w2est,w2Besti,w2Best, w2estx, w02est,w02Besti,w02Best, w02estx, extest,extBesti,extBest,extestx] =...
     marg(post,prstep,prmin,prmax,w1step,w1min,w1max,w2step,w2min,w2max,w02step,w02min,w02max,extstep,extmin,extmax); %Marginalize
    
    lw1 = length(w1est); lw2 = length(w2est); lpr = length(prest); lw02 = length(w02est); lext = length(extest);%Rename length of marginalized vectors
    sl = 4; sr = 4; % Sets how many data points to the left and right of the threshold to scan on the next loop;
    
    if l1 ==1
        [errPR_l1, errw02_l1, errw2_l1, errw1_l1] = errorcheck(prBesti,prest,prBest,w02Besti,w02est,w02Best,w2Besti,w2est, w1Besti,w1est,l1);
    end
    
    if l1==2
        [errPR_l2, errw02_l2, errw2_l2, errw1_l2] = errorcheck(prBesti,prest,prBest,w02Besti,w02est,w02Best,w2Besti,w2est, w1Besti,w1est,l1);
        
        if errPR_l1 + errw02_l1 + errw2_l1+ errw1_l1 + errPR_l2 + errw02_l2 + errw2_l2 + errw1_l2 > 0
            error = 'Yes';
        else
            error = 'No';
        end     
    end
      
    if l1 ==1
        [w1min, w1max, error1] =   param(w1est, thr, w1min, w1max, w1step, sl+20, sr+20);
        [w2min, w2max, error2] =   param(w2est, thr, w2min, w2max, w2step, sl+20, sr+20);
        [prmin, prmax, error3] =   param(prest, thr, prmin, prmax, fracstep, sl, sr);
        [w02min, w02max, error4] = param(w02est, thr, w02min, w02max, fracstep, sl, sr);
        [extmin, extmax, error5] = param(extest, thr, extmin, extmax, fracstep, sl, sr);
    end
    
    %%%%%%%%%%%%%%%% PLOTTING %%%%%%%%%%%%%
    figure((1+(l1-1)*4)); clf; plot(prestx,prest); xlabel('pr'); ylabel('P(pr)');
    figure((2+(l1-1)*4)); clf; plot(w1estx,w1est); xlabel('W1'); ylabel('P(W1)');
    figure((3+(l1-1)*4)); clf; plot(w02estx,w02est); xlabel('w02'); ylabel('P(w02)');
    figure((4+(l1-1)*4)); clf; plot(w2estx,w2est); xlabel('W2'); ylabel('P(W2)');
    
    %Plot Residue and check how good a fit we found
    if l1 ==2
        [maxpost,ind] = max(post(:));
        [w1am,w2am,w02am,pram,extam] = ind2sub(size(post),ind); %This code grabs the max values from the likelihood function
        
        w2b = (w2am-1)*w2step + w2min; %Finds the value of w2 at the max of the likilihood (different than most likely w2)
        w1b = (w1am-1)*w1step + w1min;
        w02b = (w02am - 1)*w02step + w02min;
        prb = (pram - 1)*prstep + prmin;
        extb = (extam - 1)*extstep + extmin;
        
%         f2 = exp(-s/w2b); %signal over one period
%         f2 = [f2 f2]; %signal over 2 consecutive periods
%         f2con = conv(f2,ga); %PDF after conv
%         f2bar = f2con(bins+1:2*bins); %pdf after mod-ing by laser period
%         f2h = f2bar(1:binskeep); %Keep only the appropriate bins
%         b = sum(f2h)/sum(f2bar); %Norm constant for entire period / Norm constant for window period
%         f2h = f2h/sum(f2h); %Normalized PDF of detecting a photon in a bin given photon came from w2 (non-FRET)
        
        f1 = exp(-s/w1b);
        f = [f1 f1];
        fcon = conv(f,ga);
        fbar = fcon(bins+1:2*bins);
        fh = fbar(1:binskeep);
        a = sum(fh)/sum(fbar); %Norm constant for entire period / Norm constant for window period
        fh = fh/sum(fh);
        prf = (a/b);
        
        model = wig.*((1-prb*prf*w02b-w02b-extb)/binskeep+(prb*prf*w02b)*fh + w02b*f2h + extb*ext);
        model = model/sum(model); %Normalize model
        data = p/sum(p); %Normalize data
        res = (data-model)./sqrt(data); %Caluculate Residue
        
        fprintf('\n Filename is %s irfname is %s\n',dataname, irfname);
%         fprintf(' prBest is %3.5f\n', prBest);
%         fprintf(' FRET lifetime is %3.4fns\n', w1Best);
%         fprintf(' non-FRET lifetime is %3.4fns\n\n\n', w2Best);
%         
        figure((5+(l1-1)*4)); clf; plot((data)); hold on; plot((model), 'r');
        figure((6+(l1-1)*4)); clf; plot(res);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end

%%
timestamp = datestr(clock);

xls_name = 'FLIM Analysis.xls';
[~,~,prevdata] = xlsread(xls_name);
newdata = {error, w1Best, w2Best, prBest, w02Best, [pth_data dataname], timestamp, [pth_irf irfname],...
    [pth_data_for_shift data_shift_name],w1min, w1max, w1step0, w2min, w2max,w2step0,prmin,prmax,prstep0,w02min,w02max,w2step0};
concat = [prevdata ; newdata];
xlswrite(xls_name, concat);

if strcmp(error,'Yes')
                fprintf('\n\nDANGER, Will Robinson! - You triggered an ERROR\n\n');
                [adio,Fs] = audioread('danger_WR.wav');
                sound(adio,Fs)
                pause(5);
end


% [adio,Fs] = audioread('wario.wav');
% sound(adio,Fs)

        
        prBestmat2(expt,index) = prBest;
        w2Bestmat2(expt,index) = w2Best;
        w1Bestmat2(expt,index) = w1Best;
        
    end
end

toc
 [adio,Fs] = audioread('wario.wav');
 sound(adio,Fs)
