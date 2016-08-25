%% Clears and Comments
%Here we prepare the data and set the analysis parameters

%tsm: if set to >1, appends '_cX' to file name, where X is value of tsm.
%ngr = divides FLIM data into ngr pixel groups. Ignored if 'FLIMage' is set.
%num_int_bins: Set this to the number intenisty bins you want for binning
%pixels (via binning by intensity) into SUPER PIXELS
%combine_exposures_FLIMage: number of FLIM data sets (exposures) you want
%to combine into one FLIM data.
clear;
tsm=10; %%This is for concatanating multiple images into one large image. tsm < 100;
num_int_bins = 0;
combine_exposures_FLIMage =0; %Use for
reach = 0;% Used for boxcar averaging FLIM data
ngr = 16;%jind*100000;
split_matin = 1; %Set to 1 to "split" set into one group, set to >1 for number of matins you want
save_image = 1;

set_matnum = 0;
tfw = 0;
tbac = 0;
w1vec = .25:0.05:3;% [] % Leave empty unless you want to do a w1sweep
base_name = [];
cpath = 'C:\Users\Bryan\Documents\MATLAB\data\2016-08-09\FlAsH_and_ReAsH_1to1\';
int_cor = 'alexa488_16x';

%%Cell used for the data. A new matin will be created for each filename
if ~isempty(base_name)
    %cd(cpath);
    name_list  = ls([cpath,base_name,'*']);
    %cd('C:\Users\Bryan\Documents\MATLAB\FLIM\');
    dlmwrite('datanames_preview.txt',name_list,'delimiter','');
    winopen('C:\Users\Bryan\Documents\MATLAB\FLIM\datanames_preview.txt')
    size_name_list = size(name_list);
    junk = input('Check file name list then press enter if ok');
    for dc_ind = 1:size_name_list(1)
        name=strrep(name_list(dc_ind,:),' ','');
        dataname_cell{dc_ind} = name;
    end  
end
dataname_cell ={'01','02','03','04','05'};

%dataname_cell = {'UF1_2e5_100sec','UF2_2e5_100sec','UF3_2e5_100sec','UF4_2e5_100sec_noasters'...
%                'UF1_1e6_100sec','UF2_1e6_100sec','UF3_1e6_100sec','UF4_1e6_100sec_noasters'};

%dataname_cell=  {'DR1_2e5_100sec','DR2_2e5_100sec','DR3_2e5_100sec','DR4_2e5_100sec_noasters',...
%               'DR1_1e6_100sec', 'DR2_1e6_100sec', 'DR3_1e6_100sec', 'DR4_1e6_100sec_noasters'};
                
if set_matnum
    junk = input('You are setting the matnum. Careful...Press enter to continue\n'); %#ok<*UNRCH>
end
            
%% Set search parameters
w1step = .01; w1min= .2; w1max = 3;%1.73    %1.064-2016/2/27  --- %1.73  %1.6 used for extract 8/27 and 9/5 (actual 9/5 is 1.59) %1.5 used for taxol extract; %.8-2 1.05 %w1min must be an integer multiple of w1step.
w2step = .01; w2min =  3.69; w2max = 3.69; %3.77 %3.67  %3.802     %3.82/1.49   %3.8 used for 3/7/15 data %3.745 usd for 8/27 E %3.87 used for taxol extract; %3.81 was found for 8/25 b80 and 8/24 extract

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
for dataname_ind = dataname_cell
    %[input] = ini_input(1,1,1); %Initialize input structure
    input = [];
    %for reach = 0:0
    comment = sprintf('dataname is %s',dataname_ind{1});
    for expt = 1:exptmax %determined by number of time series
        %dataname = dname(dataname,tseries, add_num_2_dataname, expt, exptmax);
        for cindex = 1:cyclesmax %determined by number of w2 spots
            %dataname = dataname;
            %% Select files for IRF, wigs, IRF shift, wigs shift, and extract signal
            dataname = dataname_ind{1};
            %cpath = 'C:\Users\Bryan\Documents\MATLAB\data\2016-07-31\';
            pth_irf = cpath; %file path IRF
            pth_data = cpath; %file path data
            pth_wigs = 'C:\Users\Bryan\Documents\MATLAB\data\2016-06-22\'; %file path wiggles
            pth_data_for_shift = cpath;
            irfname = 'irf';
            data_shift_name = 'alexa488_16x';%'tholand1_128x128_sec3'; %The IRF can be a little offset (in time) from the data, this data is used to align/find the offset and shift the data
            wigsname = 'wigs';
            pth_ext = pth_wigs; %ignore this
            extname = wigsname; %IGRNORE THIS
            
            shiftstep = .2; shiftmin = 0; shiftmax = 30;
            w2step_shift = .0025; w2min_shift = 3; w2max_shift = 4;
            backstep = .01; backmin = .01; backmax = .1;
            
            %%   Make IRF, wigs, irf shift
            sysinfo = 0; % Set to 1 if you want to force a rerun of make_irf_wig_ext
            cppout = fopen('SysInfo.txt');
            [old_filenames, count] = fscanf(cppout, '%s');
            new_filenames = [pth_irf, irfname, pth_wigs, wigsname,...
                pth_ext,extname, pth_data_for_shift, data_shift_name,...
                num2str(shiftstep),num2str(shiftmin,shiftmax), num2str(w2step_shift),...
                num2str(w2min_shift),num2str(w2max_shift),num2str(backstep),...
                num2str(backmin),num2str(backmax)];
            fclose(cppout);
            if isequal(new_filenames,old_filenames) && sysinfo == 0;
            else
                [~, ~, ~, ~, ~, ~, ~] = make_irf_wig_ext_minres_256(pth_irf, irfname, pth_wigs,...
                    wigsname, pth_ext, extname, data_shift_name, pth_data_for_shift,...
                    shiftstep, shiftmin, shiftmax, w2step_shift, w2min_shift, w2max_shift, backstep, backmin, backmax);
                fileID = fopen('SysInfo.txt','w');
                fprintf(fileID,'%s%s\n%s%s\n%s%s\n%s%s\n',pth_irf, ...
                    irfname, pth_wigs, wigsname, pth_ext,extname,pth_data_for_shift,...
                    data_shift_name,num2str(shiftstep),num2str(shiftmin,shiftmax),...
                    num2str(w2step_shift), num2str(w2min_shift),num2str(w2max_shift),...
                    num2str(backstep),num2str(backmin),num2str(backmax));
                fclose(fileID);
            end
            %% Load in IRF, wigs
            tempf=load(strcat(pth_irf,'current.mat'),'-mat','bneed',...
                'pulsewb', 'tmini', 'tmaxi', 'ext','irfname','wigsb','gab','shiftb');
            brem = tempf(1).bneed;   %Number of bins to remove from the 12.58ns to match BH data
            bins = tempf(1).pulsewb; %Number of bins that make up 12.58ns
            tmini = tempf(1).tmini;
            tmaxi = tempf(1).tmaxi;
            ext = tempf(1).ext;
            wig = tempf(1).wigsb;
            ga = tempf(1).gab; %ga is IRF after shift
            binskeep = bins-brem;
            jind =0;
            
            while jind < jmax %jind is pixel group# jmax is the total number of pixel groups
                jind = jind +1;
                %% This section is where data is saved
                %This sub-section saved FLIMages and all pixel groups.
                if jind ==1 && combine_exposures_FLIMage > 0
                    for k = 1: combine_exposures_FLIMage
                        dataname2 = strcat(dataname,'_c0',num2str(k)); %make general when cleaning up this code
                        if k==1
                            [combined_FLIMage,jmax,combined_ni] = ...
                                spc_2_his_256(tmini,tmaxi,dataname2,pth_data,...
                                ngr,tsm,'FLIMage',cpath,int_cor,reach);
                        else
                            [pmat,jmax,ni] = spc_2_his_256(tmini,tmaxi,...
                                dataname2,pth_data,ngr,tsm,'FLIMage',cpath,int_cor,reach);
                            combined_FLIMage = combined_FLIMage + pmat;
                            combined_ni=ni+combined_ni;
                        end
                    end
                    pmat = combined_FLIMage;
                    ni = combined_ni;
                    jmax=1;
                elseif jind == 1 && num_int_bins==0
                    [pmat,jmax,ni] = spc_2_his_256(tmini,tmaxi,dataname,...
                        pth_data,ngr,tsm,cpath,int_cor);
                elseif jind == 1 && num_int_bins>0
                    [pmat,~,ni] = spc_2_his_256(tmini,tmaxi,dataname,...
                        pth_data,1,tsm,'FLIMage',cpath,int_cor);
                end
                p = pmat(jind,:);
                
                %This sub-section combines pixel groups from with similar intensities into SUPER PIXELS
                %Number of photons per super pixel is saved in nph_counts. The nph_counts
                %field in the input file only exists when this section is run
                
                %Average number of photons per pixel in the super pixel is saved under nph_mean->ni.
                if jind==1 && num_int_bins>0
                    intbin = min(ni)+(0:num_int_bins)*(max(ni)-min(ni))/num_int_bins;
                    nph_mean = [];
                    dataout = [];
                    pixel_counts = [];
                    for k = 1:num_int_bins
                        nph = 0;
                        datag = 0;
                        count =0;
                        for l = 1:length(ni)
                            if intbin(k) <= ni(l) && intbin(k+1)> ni(l)
                                count = count  +1;
                                nph = ni(l)+nph;
                                datag = pmat(l,:)+datag;
                            elseif k==num_int_bins && intbin(k) <= ni(l) && intbin(k+1)>= ni(l)
                                count = count  +1;
                                nph = ni(l)+nph;
                                datag = pmat(l,:)+datag;
                            end
                        end
                        if count > 0
                            nph_mean =[nph_mean nph/count];
                            dataout = [dataout;datag];
                            pixel_counts = [pixel_counts count];
                        end
                    end
                    jmax = length(nph_mean);
                end
                
                if num_int_bins>0
                    p = dataout(jind,:);
                    ni = nph_mean;
                    input(cindex,expt,jind).pixel_counts = pixel_counts;
                    input(cindex,expt,jind).num_int_bins = num_int_bins;
                end
                
                %% Save the info into the input file %%
                for save_input_mat = 1:1
                    if split_matin==0
                        split_matin=1;
                    end
                       
                    input(cindex,expt,jind).jmax = jmax/split_matin;
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
                    
                    input(cindex,expt,jind).backstep = backstep; input(cindex,expt,jind).backmin = backmin; input(cindex,expt,jind).backmax = backmax;
                    input(cindex,expt,jind).w2step_shift = w2step_shift; input(cindex,expt,jind).w2min_shift = w2min_shift; input(cindex,expt,jind).w2max_shift = w2max_shift;
                    input(cindex,expt,jind).shiftstep = shiftstep; input(cindex,expt,jind).shiftmin = shiftmin; input(cindex,expt,jind).shiftmax = shiftmax;
                    input(cindex,expt,jind).shiftb = tempf(1).shiftb;
                    %                         input(cindex,expt,jind).r1s = r1s; input(cindex,expt,jind).r2s = r2s;
                    %                         input(cindex,expt,jind).r1l = r1l; input(cindex,expt,jind).r2l = r2l;
                    
                    input(cindex,expt,jind).dataname = dataname;
                    input(cindex,expt,jind).pth_data = pth_data;
                    input(cindex,expt,jind).irf_name= irfname;
                    input(cindex,expt,jind).pth_irf = pth_irf;
                    
                    input(cindex,expt,jind).data_shift_name = data_shift_name;
                    input(cindex,expt,jind).pth_data_for_shift = pth_data_for_shift;
                    input(cindex,expt,jind).ngr = ngr;
                    input(cindex,expt,jind).ni = ni;
                    
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
                    input(cindex,expt,jind).split_matin = split_matin;
                    input(cindex,expt,jind).reach = reach;
                    input(cindex,expt,jind).combine_exposures_FLIMage = combine_exposures_FLIMage;               
                end
            end
        end
    end
    
    if save_image
        [int_image,~,~] = spc_2_his_256(tmini,tmaxi,[dataname,'_c01'],pth_data...
            ,1,0,'save_int_image',cpath,int_cor); %#ok<UNRCH>
    else
        int_image = 0;
    end
    
    %%
    if split_matin <2        
        [set_matnum] = save_matin(input,int_image,set_matnum,dataname,w1vec);               
    else
        [set_matnum] = split_input(input,int_image,set_matnum,dataname,split_matin);
    end
end






tseries = 0; %add_num_2_dataname = 'n;
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
