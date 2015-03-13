%% Clears and Comments
%Here we prepare the data and set the analysis parameters
%pr means FRET fraction now
%clear;
loadname = 'sh-his';    
cdye = load(loadname);

% di1 = 18;
% di2 = 30;

for el = 1:1 %0:di1*di2-1
    %w2vec = [3.8201, 3.8527, 3.8627, 3.8683, 3.8727, 3.8783, 3.8730, 3.8810, 3.8910, 3.8955];
    jmax = 1;
    exptmax = 1;
    cyclesmax = 100;
    [input] = ini_input(cyclesmax,exptmax,jmax); %Initialize input structure
    comment = 'find lifetimes';
    tfw =0;
    tbac =0;
    
    for jind = 1:jmax
        for expt = 1:exptmax
            if expt == 1         
                for i=1:cyclesmax
                    str(i).name = loadname;
                end
                %                 str(2).name = 'm0n-spot2.sdt';
            end
            irfname = 'irf-347pm.sdt';
            data_shift_name = 'manual';
            
            for cindex = 1:cyclesmax
                %% Make IRF, wigs, IRF shift, wigs shift, and extract signal, read in data path
                
                sysinfo = 0; % Set to 1 if you want to force a rerun of make_irf_wig_ext
                pth_irf = 'C:\Users\Bryan\Documents\MATLAB\Bayes_2014\data\2014-10-7\'; %Type the file path of the Instrument Response Function (IRF) here.
                %irfname = 'irf-641pm.sdt';
                
                pth_wigs = 'C:\Users\Bryan\Documents\MATLAB\Bayes_2014\data\2014-10-7\';
                wigsname = 'wigs_3e5.sdt';
                
                pth_ext = 'C:\Users\Bryan\Documents\MATLAB\Bayes_2014\data\2014-6-13\';
                extname = 'extract-10m.sdt';
                
                pth_data_for_shift = pth_irf;
                
                cppout = fopen('SysInfo.txt');
                [old_filenames, count] = fscanf(cppout, '%s');
                new_filenames = [pth_irf, irfname, pth_wigs, wigsname, pth_ext,extname, pth_data_for_shift, data_shift_name];
                fclose(cppout);
                
                sfracv = 0.25; shiftmin = -10; shiftmax = 5;
                w2step_shift = 0.0025; w2min_shift = 3; w2max_shift = 4;
                backstep = 0.01; backmin = 0.001; backmax = 0.02;
                wigstep = 2; wigmin = 0; wigmax = 0;
                
                if isequal(new_filenames,old_filenames) && sysinfo == 0;
                    
                else   
                    [~, ~, ~, ~, ~, ~, ~, ~, ~] = make_irf_wig_ext(pth_irf, irfname, pth_wigs, wigsname, pth_ext, extname, data_shift_name, pth_data_for_shift, sfracv, shiftmin, shiftmax, w2step_shift,w2min_shift, w2max_shift, backstep, backmin, backmax, wigstep, wigmin, wigmax);
                    fileID = fopen('SysInfo.txt','w');
                    fprintf(fileID,'%s%s\n%s%s\n%s%s\n%s%s\n',pth_irf, irfname, pth_wigs, wigsname, pth_ext,extname, pth_data_for_shift, data_shift_name);
                    fclose(fileID);
                    %fprintf('%s', data_shift_name);
                end
                
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
                
                %% Set search parameters
                
                w1step = .005; w1min= 1; w1max = 1; %.8-2 1.05 %w1min must be an integer multiple of w1step.
                w2step = .005; w2min= cdye.info.his.w2min; w2max= cdye.info.his.w2max; %cdye 3.938% 3.745 for tfw=.6 on 8/27/14 extract %3.81
                
                fracstep = 0.002;%.0005%.0025
                prstep = fracstep; prmin=0; prmax = 0;
                w02step = fracstep; w02min = 0; w02max = 1;
                extstep = fracstep; extmin = 0; extmax = 0;
                thr = 0.01;         
                %% Set J dependence on how data is analyzed
                
                pth_data = pth_irf;
                dataname = str(cindex).name;
                thrt = 0; %0 means do not threshold, 1 means threshold , 2 means treshold for desired photon amount
                imthr = 2;%amount to threshol by. Thresh is max(int)/imthr
                DesPho = 500000;%jind*100000;
                
%                 tempdata =  load('sh-his');
%                 dtemp = tempdata.hisshort;
                %datatemp = data(el,cindex).his;
                % datatemp = data(ceil((el+1)/di2),mod(el,di2)+1).his;
                
                %p = datatemp(tmini:tmaxi)';
                p = cdye.info.his.his(tmini:tmaxi)';
                % p = data.his(tmini:tmaxi)';
                % int = 10^-6;   
                intl = cdye.info.his.int;%(1.4*10^5)/(8*10^7);%int in units of laser pulse periods             
                ints = .1;
                %  [p] = read_imagethr_14ns(tmini,tmaxi,dataname,pth_data,imthr,thrt,DesPho);
                
                %% Save the info into the input file %%              
                input(cindex,expt,jind).jmax = jmax;
                input(cindex,expt,jind).exptmax = exptmax;
                input(cindex,expt,jind).cyclesmax = cyclesmax;
                
                input(cindex,expt,jind).datahis = p;
                input(cindex,expt,jind).ga = ga; %ga is name of vector of shifted irf
                
                input(cindex,expt,jind).w1step = w1step; input(cindex,expt,jind).w1min = w1min; input(cindex,expt,jind).w1max = w1max;
                input(cindex,expt,jind).w2step = w2step; input(cindex,expt,jind).w2min = w2min; input(cindex,expt,jind).w2max = w2max;
                input(cindex,expt,jind).prstep = prstep; input(cindex,expt,jind).prmin = prmin; input(cindex,expt,jind).prmax = prmax;
                input(cindex,expt,jind).w02step = w02step; input(cindex,expt,jind).w02min = w02min; input(cindex,expt,jind).w02max = w02max;
                input(cindex,expt,jind).extstep = extstep; input(cindex,expt,jind).extmin = extmin; input(cindex,expt,jind).extmax = extmax;
                input(cindex,expt,jind).fracstep = fracstep;
                
                input(cindex,expt,jind).dataname = dataname;
                input(cindex,expt,jind).pth_data = pth_data;
                input(cindex,expt,jind) .irf_name= irfname;
                input(cindex,expt,jind).pth_irf = pth_irf;
                
                input(cindex,expt,jind).data_shift_name = data_shift_name;
                input(cindex,expt,jind).pth_data_for_shift = pth_data_for_shift;
                
                
                input(cindex,expt,jind).thrt = thrt; %0 means do not threshold data, 1 means threshold , 2 means treshold for desired photon amount
                input(cindex,expt,jind).imthr = imthr; %if input.thresh/thrt is 1, thresh is amount to threshol by. Thresh is max(int)/imthr
                input(cindex,expt,jind).DesPho = DesPho;  %if input.thresh/thrt is 2, DesPho sets the number of photons
                
                input(cindex,expt,jind).thr = thr;
                input(cindex,expt,jind).brem = brem;
                input(cindex,expt,jind).bins= bins;
                input(cindex,expt,jind).tmini = tmini;
                input(cindex,expt,jind).tmaxi = tmaxi;
                input(cindex,expt,jind).ext= ext;
                input(cindex,expt,jind).wig = wig;
                
                input(cindex,expt,jind).sfracv = sfracv; input(cindex,expt,jind).shiftmin = shiftmin; input(cindex,expt,jind).shiftmax= shiftmax;
                input(cindex,expt,jind).w2step_shift = w2step_shift; input(cindex,expt,jind).w2min_shift = w2min_shift; input(cindex,expt,jind).w2max_shift = w2max_shift;
                input(cindex,expt,jind).backstep = backstep; input(cindex,expt,jind).backmin = backmin; input(cindex,expt,jind).backmax = backmax;
                input(cindex,expt,jind).wigstep = wigstep; input(cindex,expt,jind).wigmin = wigmin; input(cindex,expt,jind).wigmax = wigmax;
                
                input(cindex,expt,jind).pth_wigs = pth_wigs;
                input(cindex,expt,jind).wigsname = wigsname;
                input(cindex,expt,jind).pth_ext = pth_ext;
                input(cindex,expt,jind).extname = extname;
                input(cindex,expt,jind).comment = comment;
                input(cindex,expt,jind).tbac = tbac;
                input(cindex,expt,jind).tfw = tfw;
                input(cindex,expt,jind).intl = intl;
                input(cindex,expt,jind).ints = ints;
                
            end
        end
    end
    [MatName,SimName] = write_to_mlist; fprintf('%s\n',MatName);
    save(MatName, 'input');
    
    if exist('sim','var') ==1
        save(SimName,'sinfo');
    end
end



%Estimate time to search space
%             sizell = size(loglike,1)*size(loglike,2)*size(loglike,3)*size(loglike,4)*size(loglike,5)/1000;
%             time_est= ceil(.06*sizell);
%fprintf('\n\n%1.0fk gridpoints. Runtime %1.0f seconds\n',sizell,time_est);
