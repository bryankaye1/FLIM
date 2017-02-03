%% Clears and Comments
%Here we prepare the data and set the analysis parameters

%tsm: if set to >1, appends '_cX' to file name, where X is value of tsm.
%ngr = divides FLIM data into ngr pixel groups. Ignored if 'FLIMage' is set.
%num_int_bins: Set this to the number intenisty bins you want for binning
%pixels (via binning by intensity) into SUPER PIXELS
%combine_exposures_FLIMage: number of FLIM data sets (exposures) you want
%to combine into one FLIM data.
clear;

set_matnum = 0;%Use this to set the matnumber
num_int_bins = 0; %Use this to create equally spaced intensity groups.
ngr = 1;%jind*100000;
split_matin = 1; %Set to 1 to "split" set into one group, set to >1 for number of matins you want

tfw = 0;
tbac = 0;
base_name = 'DONOR1_ACC9_CELL19';
dataname_cell = {'DONOR1_ACC9_CELL23'};
cpath = '/Users/bryankaye/Documents/MATLAB/data/2016-11-1/';
int_cor = 'coumarin_2e5_pdms_100sec';
data_shift_name = 'coumarin_2e5_pdms_100sec';%'DONOR_NORAN2_c99';%'uf1_2min_c50';The IRF can be a little offset (in time) from the data, this data is used to align/find the offset and shift the data
skip_remake_shift = 1;

%This section is for parameters that are zero for time-series analysis
tsm= 1; %%This is for concatanating images that all end in '_C#' into one large image. tsm < 100;
segment_FLIMdata=1; blurr = 2; im_thr = .03;

w1step = .01; w1min= 2.11; w1max = 2.11; %.97 for 11-4 extract%2.11 used for cells
w2step = .01; w2min = 3.62; w2max = 3.62; %3.62 used for cells. %3.68 used for extract

reach = 0;
make_FLIMage = 0;% Used for boxcar averaging FLIM data %Set to
combine_exposures = 0; %Used for adding exposures together
w1vec =  [];%.25:.05:2; %Set this vector to the ADDITIONAL w1 you want to set by creating new matins. Leave empty unless you want to do a w1sweep

%%Cell used for the data. A new matin will be created for each filename
if ~isempty(base_name)
[tsm,dataname_cell] = find_filenames(cpath,base_name,tsm);
end
%dataname_cell={'DONOR1_ACC9_CELL8'};
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
thr = .01; %thr is the marginalization threshold

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
            pth_wigs = cpath;%'/Users/bryankaye/Documents/MATLAB/data/2017-01-03/'; %file path wiggles
            pth_data_for_shift = cpath;
            irfname = 'from matin25255';
            wigsname = 'from matin25255';
            pth_ext = pth_wigs; %ignore this
            extname = wigsname; %IGRNORE THIS
            
            shift.step = .2; shift.min = -30; shift.max = 5; %shiftstep = .2
            shift.w2step = .025; shift.w2min = 3; shift.w2max = 4; %w2step = .0025
            shift.backstep = .01; shift.backmin = .01; shift.backmax = .1;
            
            %%   Make IRF, wigs, irf shift
            sysinfo = 0; % Set to 1 if you want to force a rerun of make_irf_wig_ext. other
            make_sysinfo(sysinfo,skip_remake_shift,cpath, pth_irf, irfname, pth_wigs,...
            wigsname, pth_ext, extname, data_shift_name, pth_data_for_shift,shift);
            
        %% Load in IRF, wigs
            load(strcat(pth_irf,'current.mat'));
            jind =0;
            while jind < jmax %jind is pixel group# jmax is the total number of pixel groups
                jind = jind +1;
                %% This section is where data is saved
                %This sub-section saved FLIMages and all pixel groups.
                if jind==1
                    if combine_exposures > 0
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
                        if tsm(1)>0 %Concatanates files ending in .c# 
                            ts = tsm(dataname_ind);
                        else
                            ts = 0;
                        end
                        [pmat,ni,int_image] = loadspc(tmini,tmaxi,dataname,...
                            pth_data,int_cor,cpath,ts,reach);
                    end
                    if segment_FLIMdata
                        input(1,1,1).blurr = blurr;
                        input(1,1,1).im_thr = im_thr;
                        [pmat,ni,segment_results] = segment_FLIMage(ni,pmat,blurr,im_thr); %normally set to 2 and .05
                    else
                        segment_results = 0;
                    end
                    if make_FLIMage
                        jmax = length(pmat);
                    elseif ngr>1
                        [pmat,jmax,ni,pixs_per_bin] = makengr(pmat,ni,ngr);
                    elseif num_int_bins>0
                        [jmax, nph_mean, pixel_counts, dataout] = make_int_bins(ni,num_int_bins,pmat);
                    elseif ngr==1
                        pmat = sum(pmat,1);
                        jmax = 1;
                        ni = sum(ni); %sum(pmat)
                    end
                end
                
                if num_int_bins>0
                    p = dataout(jind,:);
                    ni = nph_mean;
                    input(1,1,1).pixel_counts = pixel_counts;
                    input(1,1,1).num_int_bins = num_int_bins;
                else
                    p = pmat(jind,:);
                end
                
                %% Save the info into the input file %%
                [input]  = contruct_input_struct(split_matin,jmax,exptmax,...
                    cyclesmax,p,gab,w1step,w1min,w1max,w2step,w2min,w2max,prstep,prmin,prmax,...
                    w02step,w02min,w02max,fracstep,shift,shiftb,dataname,pth_data,irfname,pth_irf,...
                    data_shift_name,pth_data_for_shift,ngr,ni,thr,bneed,pulsewb,tmini,tmaxi,...
                    ext,wigsb,pth_wigs,wigsname,pth_ext,extname,comment,tbac,tfw,reach,...
                    combine_exposures,tsm,cindex,expt,jind,input);
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
