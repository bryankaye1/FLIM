%% Clears and Comments
%Here we prepare the data and set the analysis parameters
%pr means FRET fraction now

tsm=1; %%This is for concatanating multiple images into one large image. tsm < 100;
tseries = 0; findw = 3;
ngr = 1;%jind*100000;
jmax = 1;
cyclesmax = 1;
if findw > tseries
    exptmax = findw;
else
    exptmax = tseries;
end
[input] = ini_input(cyclesmax,exptmax,jmax); %Initialize input structure
comment = 'comment here';
tfw = 0.3;
tbac = 0.1;

%% Set search parameters
w1step = .02; w1min= 1.5; w1max = 1.5;%1.73  %1.6 used for extract 8/27 and 9/5 (actual 9/5 is 1.59) %1.5 used for taxol extract; %.8-2 1.05 %w1min must be an integer multiple of w1step.
w2step = .02; w2min =  2; w2max = 4;%3.802%3.82/1.49   %3.8 used for 3/7/15 data %3.745 usd for 8/27 E %3.87 used for taxol extract; %3.81 was found for 8/25 b80 and 8/24 extract

fracstep = 0.005; %.005 with w1/w2 set is 10sec per group
prstep = fracstep; prmin=0; prmax = 0;
w02step = fracstep; w02min = 0; w02max = 1;
thr = 0.01;

for expt = 1:exptmax %determined by number of time series
    dataname = 'no_acc_fov';
    dataname = dname(dataname,tseries, findw, expt, exptmax);
    for cindex = 1:cyclesmax %determined by number of w2 spots
        jind =0;
        while jind < jmax %j is determined by number of pixel groups
            jind = jind +1;
            %% Select files for IRF, wigs, IRF shift, wigs shift, and extract signal
            sysinfo = 1; % Set to 1 if you want to force a rerun of make_irf_wig_ext
            cpath = 'C:\Users\Bryan\Documents\MATLAB\data\2016-02-03\';
            pth_irf = cpath; %file path IRF
            pth_data = cpath; %file path data
            pth_wigs = 'C:\Users\Bryan\Documents\MATLAB\data\2016-02-03\'; %file path wiggles
            pth_ext = cpath; %ignore this
            pth_data_for_shift = cpath;
            
            %dataname = strcat('tholand2_128x128_sec',num2str(k));
            irfname = 'IRF';
            data_shift_name = dataname;%'tholand1_128x128_sec3'; %The IRF can be a little offset (in time) from the data, this data is used to align/find the offset and shift the data
            wigsname = 'wigs';
            extname = 'wigs'; %IGRNORE THIS
            
            %%   Make IRF, wigs, irf shift
            cppout = fopen('SysInfo.txt');
            [old_filenames, count] = fscanf(cppout, '%s');
            new_filenames = [pth_irf, irfname, pth_wigs, wigsname, pth_ext,extname, pth_data_for_shift, data_shift_name];
            fclose(cppout);
            
            if isequal(new_filenames,old_filenames) && sysinfo == 0;
            else
                [~, ~, ~, ~, ~, ~, ~] = make_irf_wig_ext_temp(pth_irf, irfname, pth_wigs,...
                    wigsname, pth_ext, extname, data_shift_name, pth_data_for_shift);
                fileID = fopen('SysInfo.txt','w');
                fprintf(fileID,'%s%s\n%s%s\n%s%s\n%s%s\n',pth_irf, irfname, pth_wigs, wigsname, pth_ext,extname, pth_data_for_shift, data_shift_name);
                fclose(fileID);
                %fprintf('%s', data_shift_name);
            end
            
            %% Load in IRF, wigs
            tempf=load(strcat(pth_irf,'current.mat'),'-mat','bneed',...
                'pulsewb', 'tmini', 'tmaxi', 'ext','irfname','wigsb','gab');
            brem = tempf(1).bneed;   %Number of bins to remove from the 12.58ns to match BH data
            bins = tempf(1).pulsewb; %Number of bins that make up 12.58ns
            tmini = tempf(1).tmini;
            tmaxi = tempf(1).tmaxi;
            ext = tempf(1).ext;
            wig = tempf(1).wigsb;
            ga = tempf(1).gab; %ga is IRF after shift
            binskeep = bins-brem;
            
            %% Set J dependence on how data is analyzed
            if jind == 1
                [pmat,jmax,ni,sinti,imagedata] = spc_2_his(tmini,tmaxi,dataname,pth_data,1,1,'imout');
                p = pmat(jind,:);
            else
                p = pmat(jind,:);
            end
            
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
            input(cindex,expt,jind).extstep = fracstep; input(cindex,expt,jind).extmin = 0; input(cindex,expt,jind).extmax = 0;
            input(cindex,expt,jind).fracstep = fracstep;
            
            input(cindex,expt,jind).dataname = dataname;
            input(cindex,expt,jind).pth_data = pth_data;
            input(cindex,expt,jind).irf_name= irfname;
            input(cindex,expt,jind).pth_irf = pth_irf;
            
            input(cindex,expt,jind).data_shift_name = data_shift_name;
            input(cindex,expt,jind).pth_data_for_shift = pth_data_for_shift;
            input(cindex,expt,jind).ngr = ngr;
            input(cindex,expt,jind).ni = ni;
            input(cindex,expt,jind).sinti = sinti;
            
            input(cindex,expt,jind).thr = thr;
            input(cindex,expt,jind).brem = brem;
            input(cindex,expt,jind).bins= bins;
            input(cindex,expt,jind).tmini = tmini;
            input(cindex,expt,jind).tmaxi = tmaxi;
            input(cindex,expt,jind).ext= ext;
            input(cindex,expt,jind).wig = wig;
            
            input(cindex,expt,jind).pth_wigs = pth_wigs;
            input(cindex,expt,jind).wigsname = wigsname;
            input(cindex,expt,jind).pth_ext = pth_ext;
            input(cindex,expt,jind).extname = extname;
            input(cindex,expt,jind).comment = comment;
            input(cindex,expt,jind).tbac = tbac;
            input(cindex,expt,jind).tfw = tfw;
            input(cindex,expt,jind).analyze = 1;
        end
    end
end
%[MatName,SimName] = write_to_mlist; fprintf('DN = %s FN = %s\n',dataname,MatName);
%save(MatName, 'input');





%% SIM DATA CODE
%% sim data
%                 if cindex ==1
%                     pr = 0;
%                 else
%                     pr = .005*2^(cindex-2);
%                 end
%                 w02 = .96;
%                 w1 = 1.5;
%                 w2 = 3.87;
%                 w03 = 1/166;
%                 w3 = .2;
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



%%%%%%%
% fprintf('DN = %s',dataname);
%     if exist('sim','var') ==1
%         save(SimName,'sinfo');
%     end


%Estimate time to search space
%             sizell = size(loglike,1)*size(loglike,2)*size(loglike,3)*size(loglike,4)*size(loglike,5)/1000;
%             time_est= ceil(.06*sizell);
%fprintf('\n\n%1.0fk gridpoints. Runtime %1.0f seconds\n',sizell,time_est);
