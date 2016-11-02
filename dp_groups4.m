%% Clears and Comments
%Here we prepare the data and set the analysis parameters

%tsm: if set to >1, appends '_cX' to file name, where X is value of tsm.
%ngr = divides FLIM data into ngr pixel groups. Ignored if 'FLIMage' is set.
%num_int_bins: Set this to the number intenisty bins you want for binning
%pixels (via binning by intensity) into SUPER PIXELS
%combine_exposures_FLIMage: number of FLIM data sets (exposures) you want
%to combine into one FLIM data.
clear;

set_matnum = 23505;%Use this to set the matnumber
num_int_bins = 16; %Use this to create equally spaced intensity groups.
ngr = 1;%jind*100000;
split_matin = 1; %Set to 1 to "split" set into one group, set to >1 for number of matins you want

tfw = 0;
tbac = 0;
base_name = 'DONOR';
cpath = 'C:\Users\Bryan\Documents\MATLAB\data\2016-11-1\';
int_cor = 'coumarin_2e5_pdms_100sec';
data_shift_name = 'DONOR1_ACC9_CELL1_c2';%'uf1_2min_c50';The IRF can be a little offset (in time) from the data, this data is used to align/find the offset and shift the data

w1step = .01; w1min= .5; w1max = .5;
w2step = .01; w2min = 3.6; w2max = 3.6; %3.6 used for cells. %3.68 used for extract
%dataname_cell = {'1DONOR_9ACC_CELL3'};

%dataname_cell = {'1DONOR_199ACC_CELL1','1DONOR_199ACC_CELL2'};

%dataname_cell = {'1DONOR_49ACC_CELL1','1DONOR_49ACC_CELL2','1DONOR_49ACC_CELL3',...
%    '1DONOR_99ACC_CELL1','1DONOR_99ACC_CELL1','1DONOR_99ACC_CELL1'};

%dataname_cell = {'1DONOR_19ACC_CELL1','1DONOR_19ACC_CELL2','1DONOR_19ACC_CELL4',...
%    '1DONOR_19ACC_CELL5_interphase'};

%dataname_cell = {'1DONOR_9ACC_CELL1','1DONOR_9ACC_CELL2','1DONOR_9ACC_CELL3',...
%    '1DONOR_9ACC_CELL4_interphase'};
%This section is for parameters that are zero for time-series analysis
tsm= 5; %%This is for concatanating images that all end in '_C#' into one large image. tsm < 100;
reach = 0;% Used for boxcar averaging FLIM data %Set to
combine_exposures = 0;
make_FLIMage = 0;
segment_FLIMdata=1;
w1vec =  [];%.25:.05:2; %Set this vector to the ADDITIONAL w1 you want to set by creating new matins. Leave empty unless you want to do a w1sweep

%%Cell used for the data. A new matin will be created for each filename
if ~isempty(base_name)
[tsm,dataname_cell] = find_filenames(cpath,base_name,tsm);
end

if set_matnum
    junk = input('You are setting the matnum. Careful...Press enter to continue\n'); %#ok<*UNRCH>
    fprintf('ok. script running:\n');
end

%% Set search parameters
fracstep = 0.001; %.005 with w1/w2 set is 10sec per group
if w2min~=w2max
    prstep = fracstep; prmin=0; prmax = 0;
else
    prstep = fracstep; prmin=0; prmax = 1;
end
w02step = fracstep; w02min = 0; w02max = 1;
thr = 0.01;


jmax = 1;
cyclesmax = 1;
exptmax= 1;
dataname_matnum = {};
for dataname_ind = 1:length(dataname_cell)
    input = [];
    comment = sprintf('dataname is %s',dataname_cell{dataname_ind});
    for expt = 1:exptmax %determined by number of time series
        for cindex = 1:cyclesmax %determined by number of w2 spots
            %% Select files for IRF, wigs, IRF shift, wigs shift, and extract signal
            dataname = dataname_cell{dataname_ind};
            %cpath = 'C:\Users\Bryan\Documents\MATLAB\data\2016-07-31\';
            pth_irf = cpath; %file path IRF
            pth_data = cpath; %file path data
            pth_wigs = 'C:\Users\Bryan\Documents\MATLAB\data\2016-06-22\'; %file path wiggles
            pth_data_for_shift = cpath;
            irfname = 'irf';
            wigsname = 'wigs';
            pth_ext = pth_wigs; %ignore this
            extname = wigsname; %IGRNORE THIS
            
            shiftstep = .2; shiftmin = -15; shiftmax = 15; %shiftstep = .2
            w2step_shift = .0025; w2min_shift = 3; w2max_shift = 4; %w2step = .0025
            backstep = .01; backmin = .01; backmax = .1;
            %%   Make IRF, wigs, irf shift
            sysinfo = 0; % Set to 1 if you want to force a rerun of make_irf_wig_ext. other
            if exist([cpath,'SysInfo.txt'],'file')==2
                cppout = fopen([cpath,'SysInfo.txt'],'r+');
                [old_filenames, count] = fscanf(cppout, '%s');
                fclose(cppout);
            else
                old_filenames = 'sysinfo did not exist';
            end
            new_filenames = [pth_irf, irfname, pth_wigs, wigsname,...
                pth_ext,extname, pth_data_for_shift, data_shift_name,...
                num2str(shiftstep),num2str(shiftmin),num2str(shiftmax), num2str(w2step_shift),...
                num2str(w2min_shift),num2str(w2max_shift),num2str(backstep),...
                num2str(backmin),num2str(backmax)];
            if ~isequal(new_filenames,old_filenames) || sysinfo ~= 0;
                [~, ~, ~, ~, ~, ~, ~] = make_irf_wig_ext_minres(pth_irf, irfname, pth_wigs,...
                    wigsname, pth_ext, extname, data_shift_name, pth_data_for_shift,...
                    shiftstep, shiftmin, shiftmax, w2step_shift, w2min_shift, w2max_shift, backstep, backmin, backmax);
                fileID = fopen([cpath,'SysInfo.txt'],'w');
                fprintf(fileID,'%s%s\n%s%s\n%s%s\n%s%s\n',pth_irf, ...
                    irfname, pth_wigs, wigsname, pth_ext,extname,pth_data_for_shift,...
                    data_shift_name,num2str(shiftstep),num2str(shiftmin),num2str(shiftmax),...
                    num2str(w2step_shift), num2str(w2min_shift),num2str(w2max_shift),...
                    num2str(backstep),num2str(backmin),num2str(backmax));
                fclose(fileID);
            end
            %% Load in IRF, wigs
            load(strcat(pth_irf,'current.mat'));
            jind =0;
            while jind < jmax %jind is pixel group# jmax is the total number of pixel groups
                jind = jind +1;
                %% This section is where data is saved
                %This sub-section saved FLIMages and all pixel groups.
                if jind==1
                    if combine_exposures > 0
                        %combine_alarm= input('STOP!!! combine FLIMage not setup yet\n');
                        pmat = 0;
                        int_image = 0;
                        ni = 0;
                        for ce_ind = 1:combine_exposures
                            dataname_comexp = append_timeseries_names(dataname,ce_ind,combine_exposures);
                            [pmat_temp,ni_temp,int_image_temp] = loadspc(tmini,tmaxi,dataname_comexp,...
                                pth_data,int_cor,cpath,0,reach);
                            pmat = pmat + pmat_temp;
                            int_image = int_image + int_image_temp;
                            ni = ni + ni_temp;
                        end
                    else
                        if tsm(1)>0
                            ts = tsm(dataname_ind);
                        else
                            ts = 0;
                        end
                        [pmat,ni,int_image] = loadspc(tmini,tmaxi,dataname,...
                            pth_data,int_cor,cpath,ts,reach);
                    end
                    if segment_FLIMdata
                        [pmat,ni,segment_results] = segment_FLIMage(ni,pmat,2,.05); %normally set to 2 and .05
                    else
                        segment_results = 0;
                    end
                    if reach>0 || make_FLIMage
                    elseif ngr>1
                        [pmat,jmax,ni,pixs_per_bin] = makengr(pmat,ni,ngr);
                    elseif num_int_bins>0
                        [jmax, nph_mean, pixel_counts, dataout] = make_int_bins(ni,num_int_bins,pmat);
                    elseif ngr==1
                        pmat = sum(pmat,1);
                        jmax = 1;
                        ni = sum(pmat);
                    end
                end
                
                if num_int_bins>0
                    p = dataout(jind,:);
                    ni = nph_mean;
                    input(cindex,expt,jind).pixel_counts = pixel_counts;
                    input(cindex,expt,jind).num_int_bins = num_int_bins;
                else
                    p = pmat(jind,:);
                end
                
                %% Save the info into the input file %%
                for save_input_mat = 1:1
                    if split_matin==0
                        split_matin=1;
                    end
                    
                    input(1,1,1).jmax = jmax/split_matin;
                    input(1,1,1).exptmax = exptmax;
                    input(1,1,1).cyclesmax = cyclesmax;
                    
                    input(cindex,expt,jind).datahis = p;
                    input(cindex,expt,jind).ga = gab; %ga is name of vector of shifted irf
                    
                    input(cindex,expt,jind).w1step = w1step; input(cindex,expt,jind).w1min = w1min; input(cindex,expt,jind).w1max = w1max;
                    input(cindex,expt,jind).w2step = w2step; input(cindex,expt,jind).w2min = w2min; input(cindex,expt,jind).w2max = w2max;
                    input(cindex,expt,jind).prstep = prstep; input(cindex,expt,jind).prmin = prmin; input(cindex,expt,jind).prmax = prmax;
                    input(cindex,expt,jind).w02step = w02step; input(cindex,expt,jind).w02min = w02min; input(cindex,expt,jind).w02max = w02max;
                    input(cindex,expt,jind).extstep = fracstep; input(cindex,expt,jind).extmin = 0; input(cindex,expt,jind).extmax = 0;
                    input(cindex,expt,jind).fracstep = fracstep;
                    
                    input(1,1,1).backstep = backstep; input(1,1,1).backmin = backmin; input(1,1,1).backmax = backmax;
                    input(1,1,1).w2step_shift = w2step_shift; input(1,1,1).w2min_shift = w2min_shift; input(1,1,1).w2max_shift = w2max_shift;
                    input(1,1,1).shiftstep = shiftstep; input(1,1,1).shiftmin = shiftmin; input(1,1,1).shiftmax = shiftmax;
                    input(1,1,1).shiftb = shiftb;
                    %                         input(cindex,expt,jind).r1s = r1s; input(cindex,expt,jind).r2s = r2s;
                    %                         input(cindex,expt,jind).r1l = r1l; input(cindex,expt,jind).r2l = r2l;
                    input(cindex,expt,jind).dataname = dataname;
                    input(cindex,expt,jind).pth_data = pth_data;
                    input(cindex,expt,jind).irf_name= irfname;
                    input(cindex,expt,jind).pth_irf = pth_irf;
                    
                    input(1,1,1).data_shift_name = data_shift_name;
                    input(1,1,1).pth_data_for_shift = pth_data_for_shift;
                    input(1,1,1).ngr = ngr;
                    input(1,1,1).ni = ni;
                    
                    input(1,1,1).thr = thr;
                    input(1,1,1).brem = bneed;
                    input(1,1,1).bins= pulsewb;
                    input(1,1,1).tmini = tmini;
                    input(1,1,1).tmaxi = tmaxi;
                    input(1,1,1).ext= ext;
                    input(1,1,1).wig = wigsb;
                    
                    input(1,1,1).pth_wigs = pth_wigs;
                    input(1,1,1).wigsname = wigsname;
                    input(1,1,1).pth_ext = pth_ext;
                    input(1,1,1).extname = extname;
                    input(1,1,1).comment = comment;
                    input(1,1,1).tbac = tbac;
                    input(1,1,1).tfw = tfw;
                    input(1,1,1).split_matin = split_matin;
                    input(1,1,1).reach = reach;
                    input(1,1,1).combine_exposures_FLIMage = combine_exposures;
                    input(1,1,1).tsm = tsm;
                end
            end
        end
    end
    %%
    if split_matin <2
        [set_matnum] = save_matin(input,int_image,segment_results,set_matnum,dataname,w1vec);
        dataname_matnum{end+1} = {dataname, set_matnum-1};
    else
        [set_matnum] = split_input(input,int_image,segment_results,set_matnum,dataname,split_matin);
    end
end
%%
[start_nums,end_nums,tseries_names] = find_series(dataname_matnum);
if ~isempty(start_nums)
    for i = 1:length(tseries_names)
        fprintf( '%s_1-%s.  %3.0f:%3.0f\n', tseries_names{i}, num2str(end_nums(i)-start_nums(i)+1),...
            start_nums(i), end_nums(i));
    end
    
    fprintf('\n\n vector to enter into FFcomp script:\n [');
    for i = 1:length(tseries_names)-1
        fprintf( '%3.0f:%3.0f,', start_nums(i), end_nums(i));
    end
    fprintf( '%3.0f:%3.0f', start_nums(end), end_nums(end));
    fprintf('];\n');
end

% dataname_cell ={'ACC_ONLY_CELL1','ACC_ONLY_CELL2',...
%     '1DONOR_3ACC_CELL1','1DONOR_3ACC_CELL2','1DONOR_3ACC_CELL3_interphase',...
%    '1DONOR_1ACC_CELL1','1DONOR_1ACC_CELL2','1DONOR_1ACC_CELL3_interphase',...
%    '3DONOR_1ACC_CELL1','3DONOR_1ACC_CELL2_no_labeled_tub','3DONOR_1ACC_CELL3','3DONOR_1ACC_CELL4_interphase'...
%    'DONOR_ONLY_CELL1','DONOR_ONLY_CELL2','DONOR_ONLY_CELL3','DONOR_ONLY_CELL4','DONOR_ONLY_CELL5_interphase'};



%w1step = .01; w1min= 3; w1max = 3;%1.73    %1.064-2016/2/27  --- %1.73  %1.6 used for extract 8/27 and 9/5 (actual 9/5 is 1.59) %1.5 used for taxol extract; %.8-2 1.05 %w1min must be an integer multiple of w1step.
%w2step = .01; w2min = 1; w2max = 5;%3.69 %3.77 %3.67  %3.802     %3.82/1.49   %3.8 used for 3/7/15 data %3.745 usd for 8/27 E %3.87 used for taxol extract; %3.81 was found for 8/25 b80 and 8/24 extract




%tseries = 0; %add_num_2_dataname = 'n;
% if add_num_2_dataname > tseries
%     exptmax = add_num_2_dataname;
% else
%     exptmax = tseries;
% end
%[input] = ini_input(cyclesmax,exptmax,jmax); %Initialize input structure


%         [MatName,SimName] = write_to_mlist(set_matnum);
%         fprintf('\nDN = %s FN = %s\n',dataname,MatName);
%         fileID = fopen('matin_prints.txt','at');
%         fprintf(fileID,'DN = %s FN = %s\n',dataname,MatName);
%         fclose(fileID);
%         save(MatName, 'input','int_image');
%         if ~isempty(w1vec)
%             set_matnum = change_lifetimes(MatName,w1vec,set_matnum,int_image);
%         elseif set_matnum
%             set_matnum = set_matnum + 1;
%         end






%         input_holdon = input;
%         for k = 1:split_matin
%             pstart = 1+(k-1)*(jmax/split_matin);
%             pend = k*(jmax/split_matin);
%             input = input_holdon(1,1,round(pstart):round(pend));
%             [MatName,SimName] = write_to_mlist(set_matnum);
%             set_matnum = set_matnum + 1;
%             if k==1 % This section prints out the matin #s
%                 matstart = MatName(35:strfind(MatName,'.')-1);
%                 matend = num2str(str2num(matstart) + split_matin - 1);
%                 fprintf('DN = %s split_matin: Matin %s-%s\n',dataname,matstart,matend);
%                 fileID = fopen('matin_prints.txt','at');
%                 fprintf(fileID,'DN = %s split_matin: Matin %s-%s\n',dataname,matstart,matend);
%                 fclose(fileID);
%             end
%             save(MatName, 'input','int_image');
%         end










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