%% Clears and Comments
%Highly simplified code adapted from dp_groups3 by Bryan Kaye
clear;
ngr =1; tsm=1;
w2step = .001; w2min = .01; w2max = .5;
backstep = .001; backmin =0; backmax = 0; % For your short lifetime dye, we can 
%estimate the background by looking at the end of the FLIM decay curve


%% Select files paths...
cpath = 'C:\Users\Bryan\Documents\MATLAB\data\2015-08-27\';
pth_irf = cpath; %file path IRF
pth_data = cpath; %file path data
pth_wigs = cpath; %file path wiggles
pth_ext = cpath; %ignore this
pth_data_for_shift = cpath;

dataname = 'tholand2_128x128_fixedpoint-midplane';
irfname = 'IRF';
data_shift_name = 'tholand2_128x128_fixedpoint-midplane'; %The IRF can be a little offset (in time) from the data, this data is used to align find the offset and shift the data
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
[pmat,jmax,ni,sinti] = spc_2_his(tmini,tmaxi,dataname,pth_data,ngr,tsm);
p = pmat(1,:);

%%
Tlaser=12.58;
sumres = inf;
            
for w2i = w2min:w2step:w2max
    s = Tlaser/bins:Tlaser/bins:Tlaser;
    f2 = exp(-s/w2i); %signal over one period
    f2 = [f2 f2]; %signal over 2 consecutive periods
    f2con = conv(f2,ga); %PDF after conv
    f2bar = f2con(bins+1:2*bins); %pdf after mod-ing by laser period
    f2h = f2bar(1:binskeep); %Keep only the appropriate bins
    f2h = f2h/sum(f2h);
    for backi = backmin:backstep:backmax
        %back = backi/bins;
        back = mean(p(2500:end))/sum(p); %Here we estimate the background from FLIM decay curve
        model = (f2h + back).*wig;    
        model = model/sum(model);
        data = p/sum(p);
        res = (data-model)./sqrt(model);
        sumresi = sum(res.^2);       
        if sumresi < sumres
            sumres = sumresi;
            w2b = w2i;
            backb = backi;   
            modelb = model;
            resb = res;
        end
    end
end    
    %% Plots & irf
    
    fprintf('lifetime: %3.3f  background fraction: %3.6f Kres %3.6f \n', w2b, 100*backb, 1000*sum(res.^2));
    figure(1); clf; plot((data)); hold on; plot((modelb), 'r'); title('Data (blue) and single exponential fit (red)');
    figure(2); clf; plot(resb); title('Residue');
    figure(1); clf; hold on; 
    plot((log(data))); plot((log(modelb)), 'r'); plot(log(ga),'g');
    title('LOG: Data (blue) and single exponential fit (red)');






