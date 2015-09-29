%% Clears and Comments
%Highly simplified code adapted from dp_groups3 by Bryan Kaye
clear;

%%Set search space parameters here
w2step = .01; w2min = .1; w2max = .4;
w1step = .01; w1min = .01; w1max = .1;
astep = .01; amin = 0; amax =1;
backstep = .001; backmin =0; backmax = 0; % For your short lifetime dye, we can
%estimate the background by looking at the end of the FLIM decay curve
%%
set(0, 'DefaulttextInterpreter', 'none');
for k = 1:3
    %% Select files paths...
    cpath = 'C:\Users\Bryan\Documents\MATLAB\data\2015-08-27\';
    pth_irf = cpath; %file path IRF
    pth_data = cpath; %file path data
    pth_wigs = cpath; %file path wiggles
    pth_ext = cpath; %ignore this
    pth_data_for_shift = cpath;
    
    dataname = strcat('tholand2_128x128_sec',num2str(k));
    irfname = 'IRF';
    data_shift_name = dataname;%'tholand1_128x128_sec3'; %The IRF can be a little offset (in time) from the data, this data is used to align/find the offset and shift the data
    wigsname = 'wigs';
    extname = 'wigs'; %IGRNORE THIS
    
    %%   Make IRF, wigs, irf shift
    sysinfo=0; %Set to 1 if you want to force a rerun of IRF-shift
    cppout = fopen('SysInfo.txt');
    [old_filenames, count] = fscanf(cppout, '%s');
    new_filenames = [pth_irf, irfname, pth_wigs, wigsname, pth_ext,extname, pth_data_for_shift, data_shift_name];
    fclose(cppout);
    
    if isequal(new_filenames,old_filenames) && sysinfo == 0;
    else
        [~, ~, ~, ~, ~, ~, ~, ~, ~] = make_irf_wig_ext_temp(pth_irf, irfname, pth_wigs,...
            wigsname, pth_ext, extname, data_shift_name, pth_data_for_shift);
        fileID = fopen('SysInfo.txt','w');
        fprintf(fileID,'%s%s\n%s%s\n%s%s\n%s%s\n',pth_irf, irfname, pth_wigs, wigsname, pth_ext,extname, pth_data_for_shift, data_shift_name);
        fclose(fileID);
        %fprintf('%s', data_shift_name);
    end
    
    %% Load in IRF, wigs
    tempf=load(strcat(pth_irf,'current.mat'),'-mat','bneed',...
        'pulsewb', 'irf', 'irfsim', 'irf_pdf', 'tmini', 'tmaxi', 'ext','irfname','wigsb','gab');
    brem = tempf(1).bneed;   %Number of bins to remove from the 12.5ns to match BH data
    bins = tempf(1).pulsewb; %Number of bins that make up 12.5ns
    tmini = tempf(1).tmini;
    tmaxi = tempf(1).tmaxi;
    ext = tempf(1).ext;
    wig = tempf(1).wigsb';
    ga = tempf(1).gab/sum(tempf(1).gab);
    binskeep = bins-brem;
    
    %% Load data
    [pmat,jmax,ni,sinti, intim] = spc_2_his(tmini,tmaxi,dataname,pth_data,1,1,'imout');
    p = pmat(1,:);
    
    %%
    Tlaser=12.58;
    sumres = inf;
    s = Tlaser/bins:Tlaser/bins:Tlaser;
    for w2i = w2min:w2step:w2max
        f2t = exp(-s/w2i); %signal over one period
        f2 = [f2t f2t]; %signal over 2 consecutive periods
        f2con = conv(f2,ga); %PDF after conv
        f2bar = f2con(bins+1:2*bins); %pdf after mod-ing by laser period
        f2h = f2bar(1:binskeep); %Keep only the appropriat14e bins
        f2h = f2h/sum(f2h);
        for w1i = w1min:w1step:w1max
            f1t = exp(-s/w1i); %signal over one period
            f1 = [f1t f1t]; %signal over 2 consecutive periods
            f1con = conv(f1,ga); %PDF after conv
            f1bar = f1con(bins+1:2*bins); %pdf after mod-ing by laser period
            f1h = f1bar(1:binskeep); %Keep only the appropriate bins
            f1h = f1h/sum(f1h);
            for a = amin:astep:amax
                fs = a*f1h+(1-a)*f2h;
                fs = fs/sum(fs);
                for backi = backmin:backstep:backmax
                    %back = backi/bins;
                    back = mean(p(2500:end))/sum(p); %Here we estimate the background from FLIM decay curve
                    model = (fs + back).*wig;
                    model = model/sum(model);
                    data = p/sum(p);
                    res = (data-model)./sqrt(model);
                    sumresi = sum(res.^2);
                    if sumresi < sumres
                        sumres = sumresi;
                        w2b = w2i;
                        w1b = w1i;
                        ab = a;
                        backb = backi;
                        modelb = model;
                        resb = res;
                    end
                end
            end
        end
    end

    %% Plots & irf
    
    fprintf('short/long lifetime: %3.3f/%3.3f ns Relative Amplitude %3.2f \n\n', w1b, w2b, ab);
    figure(1); clf; plot((p)); hold on; plot((modelb*sum(p)), 'r'); title('Data (blue) and Fit (red)');
    annotation('textbox', [0.5,0.6,0.1,0.1],'String',sprintf('short/long lifetime: %3.3f/%3.3f ns',w2b,w1b));
    figure(2); clf; plot(resb); title('Residue');
    figure(3); clf; hold on;
    plot((log(p))); plot((log(modelb*sum(p))), 'r');
    annotation('textbox', [0.5,0.6,0.1,0.1],'String',sprintf('short/long lifetime: %3.3f/%3.3f ns',w2b,w1b));
    title('LOG: Data (blue) and Fit (red)');% plot(log(ga),'g');
    figure(4); clf; imagesc(squeeze(intim)); colormap('gray');
    title(dataname);
    drawnow
end






