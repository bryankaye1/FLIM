%% Clears and Comments
%Here we prepare the data and set the analysis parameters
%pr means FRET fraction now
clear;


tic
for el = 1:128 %128
    %w2vec = [3.8201, 3.8527, 3.8627, 3.8683, 3.8727, 3.8783, 3.8730, 3.8810, 3.8910, 3.8955];
    jmax = 1;
    exptmax = 1;
    cyclesmax = 128;
    [input] = ini_input(cyclesmax,exptmax,jmax); %Initilize input structure
    comment = 'changes in ff with w3 as a function of tfw';
    tfw = 0;
    tbac = 0;

    for jind = 1:jmax
        for expt = 1:exptmax
            for cindex = 1:cyclesmax
                if cindex ==1
                    
                    str(1).name = strcat('m3-spot1','.sdt');                    
                    irfname = 'irf-642pm.sdt';
                    data_shift_name = 'm0-spot1.sdt';           
                    
                    %% Make IRF, wigs, IRF shift, wigs shift, and extract signal, read in data path
                    
                    sysinfo = 0; % Set to 1 if you want to force a rerun of make_irf_wig_ext
                    pth_irf = 'C:\Users\Bryan\Documents\MATLAB\Bayes_2014\data\2014-6-13\'; %Type the file path of the Instrument Response Function (IRF) here.                    
                    pth_wigs = 'Z:\bkaye\Bayes_2013\data\2014-2-13\';
                    wigsname = 'wigs_1p5e6.sdt';                   
                    pth_ext = 'C:\Users\Bryan\Documents\MATLAB\Bayes_2014\data\2014-6-13\';
                    extname = 'extract-10m.sdt';          
                    pth_data_for_shift = pth_irf;
                    
                    cppout = fopen('SysInfo.txt');
                    [old_filenames, count] = fscanf(cppout, '%s');
                    new_filenames = [pth_irf, irfname, pth_wigs, wigsname, pth_ext,extname, pth_data_for_shift, data_shift_name];
                    fclose(cppout);
                    
                    sfracv = 0.5; shiftmin = -10; shiftmax = 5;
                    w2step_shift = 0.01; w2min_shift = 3.4; w2max_shift = 3.9;
                    backstep = 0.01; backmin = 0.01; backmax = 0.1;
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
                    
                    w1step = .005; w1min= 1.5; w1max = 1.5; %.8-2 1.05 %w1min must be an integer multiple of w1step.
                    w2step = .005; w2min =  3.87; w2max =  3.87; %3.81
                    
                    fracstep = 0.001;%.002
                    prstep = fracstep; prmin=0; prmax = 1;
                    w02step = fracstep; w02min = 0; w02max = 1;
                    extstep = fracstep; extmin = 0; extmax = 0;
                    thr = 0.01;
                    
                    
                    %% Set J dependence on how data is analyszed
                    
                    pth_data = pth_irf;
                    dataname = str(1).name;
                    thrt = 0; %0 means do not threshold, 1 means threshold , 2 means treshold for desired photon amount
                    imthr = 2;%amount to threshol by. Thresh is max(int)/imthr
                    DesPho = 500000;%jind*100000;
                    
                          
                input(1,1,1).jmax = jmax;
                input(1,1,1).exptmax = exptmax;
                input(1,1,1).cyclesmax = cyclesmax;
                input(1,1,1).ga = ga; %ga is name of vector of shifted irf                
                input(1,1,1).w1step = w1step; input(1,1,1).w1min = w1min; input(1,1,1).w1max = w1max;
                input(1,1,1).w2step = w2step; input(1,1,1).w2min = w2min; input(1,1,1).w2max = w2max;
                input(1,1,1).prstep = prstep; input(1,1,1).prmin = prmin; input(1,1,1).prmax = prmax;
                input(1,1,1).w02step = w02step; input(1,1,1).w02min = w02min; input(1,1,1).w02max = w02max;
                input(1,1,1).extstep = extstep; input(1,1,1).extmin = extmin; input(1,1,1).extmax = extmax;
                input(1,1,1).fracstep = fracstep;                
                input(1,1,1).dataname = dataname;
                input(1,1,1).pth_data = pth_data;
                input(1,1,1) .irf_name= irfname;
                input(1,1,1).pth_irf = pth_irf;                
                input(1,1,1).data_shift_name = data_shift_name;
                input(1,1,1).pth_data_for_shift = pth_data_for_shift;                          
                input(1,1,1).thrt = thrt; %0 means do not threshold data, 1 means threshold , 2 means treshold for desired photon amount
                input(1,1,1).imthr = imthr; %if input.thresh/thrt is 1, thresh is amount to threshol by. Thresh is max(int)/imthr
                input(1,1,1).DesPho = DesPho;  %if input.thresh/thrt is 2, DesPho sets the number of photons               
                input(1,1,1).thr = thr;
                input(1,1,1).brem = brem;
                input(1,1,1).bins= bins;
                input(1,1,1).tmini = tmini;
                input(1,1,1).tmaxi = tmaxi;
                input(1,1,1).ext= ext;
                input(1,1,1).wig = wig;                
                input(1,1,1).sfracv = sfracv; input(1,1,1).shiftmin = shiftmin; input(1,1,1).shiftmax= shiftmax;
                input(1,1,1).w2step_shift = w2step_shift; input(1,1,1).w2min_shift = w2min_shift; input(1,1,1).w2max_shift = w2max_shift;
                input(1,1,1).backstep = backstep; input(1,1,1).backmin = backmin; input(1,1,1).backmax = backmax;
                input(1,1,1).wigstep = wigstep; input(1,1,1).wigmin = wigmin; input(1,1,1).wigmax = wigmax;            
                input(1,1,1).pth_wigs = pth_wigs;
                input(1,1,1).wigsname = wigsname;
                input(1,1,1).pth_ext = pth_ext;
                input(1,1,1).extname = extname;
                input(1,1,1).comment = comment;
                input(1,1,1).tbac = tbac;
                input(1,1,1).tfw = tfw;   

                    
                    
                end
                if cindex ==1 && el==1
                    [p] = read_imagethr_14ns_temp(tmini,tmaxi,dataname,pth_data,imthr,thrt,DesPho);
                end
                %% sim data
                %
                %                 if cindex ==1
                %                     pr = 0;
                %                 else
                %                 pr = .005*2^(cindex-2);
                %                 end
                %                 w02 = .96;
                %                 w1 = 1.5;
                %                 w2 = 3.87;
                %
                %                 w03 = 1/166;
                %                 w3 = .2;
                %
                %                 nps = 10^7;
                %                 w01 = pr*w02*w1/(w2*(1-pr));
                %                 w00 = 1 - w01 - w02 - w03; %pr is now fret fraction
                %                 %%
                %                 [w00out, w01out, w02out, npho, p] = SimData_v2(w03,w3,pr, w02, w1, w2, nps);
                %                 sinfo(cindex,expt,jind).pr = pr;
                %                 sinfo(cindex,expt,jind).w02 = w02;
                %                 sinfo(cindex,expt,jind).w1 = w1;
                %                 sinfo(cindex,expt,jind).w2= w2;
                %                 sinfo(cindex,expt,jind).w03 = w03;
                %                 sinfo(cindex,expt,jind).w3= w3;
                %                 sinfo(cindex,expt,jind).nps = nps;
                
                %% Save the info into the input file %%
                                input(cindex,expt,jind).datahis = p(:,el,cindex)';               
            end
        end
    end
    if el==1
    [MatName,SimName] = write_to_mlisti;
    end
    nMatName = strcat(MatName(1:end-4),'-',num2str(el),MatName(end-3:end));
    save(nMatName, 'input');
    fprintf('%s\n',nMatName);
    
    if exist('sim','var') ==1
        save(SimName,'sinfo');
    end
end

toc

%Estimate time to search space
%             sizell = size(loglike,1)*size(loglike,2)*size(loglike,3)*size(loglike,4)*size(loglike,5)/1000;
%             time_est= ceil(.06*sizell);
%fprintf('\n\n%1.0fk gridpoints. Runtime %1.0f seconds\n',sizell,time_est);
