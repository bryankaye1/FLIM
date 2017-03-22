%[bneed, pulsewb, irf, irfsim, irf_pdf, wig, tmini, tmaxi, ext,varargout] =...
%    make_irf_wig_ext_minres(pth_irf, irfname, pth_wigs, wigsname, pth_ext, extname, data_shift_name, pth_data_for_shift,...
%    sfrc, shiftmin, shiftmax, w2step, w2min, w2max, backstep, backmin, backmax)

%This file is used to make the IRF vector, wiggles vector, extract vector. It also calculates the
%amount of missing time. It also makes the IRF_pdf (which is used for
%simulating data).

clear;
lifetimes = 2;
set_w00i = .02; %input 'NA" or other char array to not set w00i
dataname = 'M5';
wigsname = 'wigs_2017-03-21';
cpath = '/Users/bryankaye/Documents/MATLAB/data/2017-03-16/';

pth_wigs = cpath;
pth_data = cpath;

load(strcat(pth_data,'current.mat'));
bins = pulsewb;
Tlaser = 14;
sfrc = .2;

%% IRF and Wig shift
%Checks if files passed is in time series, if so, uses loads whole series
tsm = length(dir([pth_data,dataname,'_c*']));
[data1,~,~] = spc_2_his(tmini,tmaxi,dataname,pth_data,1,tsm);
data = data1/sum(data1);

[wigsb_un,~,~] = spc_2_his(tmini,tmaxi,wigsname,pth_wigs,1,1);
wig1 = wigsb_un/mean(wigsb_un);
m = 8;
wig =  tsmovavg(wig1,'s',m)';
wigsb = [wig1(1:m-1), wig(m:end)'];
%wigsb =  [wig; ones(bneed,1)];  %add bins (with 1 in each added bin) to make up 1 period

figure(4); plot(data);
figure(5); plot(data./wigsb);

if lifetimes==1
w2i = 3.5;
w02i = .95;
x0 =[w2i w02i];

lb = [0,0];
ub = [5,1];

elseif lifetimes==2

w2i = 3.5;
w02i = .95;
w1i = 1;
if ischar(set_w00i)
    w01i = 1-w02i;
else
    w01i = 1-w02i-set_w00i;
end
x0 =[w2i w02i w1i w01i];

lb = [.1,0,.1,0];
ub = [5,1,5,1];
end

fn_handle = @(x)lifetime_fit(x, bins, bneed, irf, Tlaser, wigsb, data,...
    lifetimes,set_w00i);
%options = optimoptions(@simulannealbnd,'FunctionTolerance',1e-50);
options = optimoptions(@simulannealbnd, 'ReannealInterval',10,'MaxTime',300,...
                     'PlotFcn',{@saplotbestf,@saplottemperature,@saplotf,...
                      @saplotstopping});    
                 
                    %'AnnealingFcn',@annealingboltz);

[x,fval,exitFlag,output] = simulannealbnd(fn_handle,x0,lb,ub,options);

fprintf('The number of iterations was : %d\n', output.iterations);
fprintf('The number of function evaluations was : %d\n', output.funccount);
fprintf('The best function value found was : %g\n', fval);
%%
[fval,model,dataplot] = fn_handle(x);
figure(2); clf; hold on;
plot(dataplot,'b');
plot(model,'r');


    
    %% Plots & irf 
%     fprintf('shift %3.4f  lifetime %3.3f  back %3.6f Kres %3.6f \n', shiftb, w2b, backb, 1000*sum(res.^2));
%     figure(11); clf; 
%     subplot(2,1,1); hold on; plot(data); plot(modelb, 'r'); title('Shift: Model vs Data');
%     subplot(2,1,2); plot(resb); title('Residue'); drawnow;


