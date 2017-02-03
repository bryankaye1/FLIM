%%Used to plot polymer series plots after reviewers comments:
% Error in the FLIM measurements of polymer is calculated from
% the standard deviation of each series measurements. 
% The normalization constant (converting num of total photons to number of
% labeled donors) is given by the highest polymer sample taken from
% matout 3229-3258. See line 167.


clear;
figure(1);clf;
imin = [3139,3154,3169,3199,3229];%
imax = [3153,3168,3198,3228,3258]; %
for loop = 1:length(imax)
    ind = 0; %index reset to zero for each dilution series
    %%Load in data for one dilution series, fit FF v int,
    %%and calculate polymer amount given fit parameters
    for i = imin(loop):imax(loop)  %Load in matout files
        clear output intb ffP cintpr
        [~,output] = load_mat_data(i);
        ind =ind+1;
        %%Set fit paramters for FF v int over image
        ni = output(1,1,1).ni; 
        al = output(1,1,1).w1Best./output(1,1,1).w2Best;
        ep = output(1,1,1).sinti; 
        %%Get FRET fraction (fluoro population) and intensity
        for j = 1:length(output)
            x(j) = (ni(j)/ep(j))*(sum(output(1,1,j).prest'.*output(1,1,j).prestx)...
               +sum(output(1,1,j).w02est'.*output(1,1,j).w02estx));
            [y(j),stdpr(j)] = transform_wf_to_f(output(1,1,j).prest,...
                output(1,1,j).prestx,al,output(1,1,j).w02est,...
                output(1,1,j).w02estx,x(j),'dont_combine');
        end
        %Remove pixel groups with zero photons
        if min(x)==0
            istart = find(x==0,1,'last')+1;
            x = x(istart:end);
            y = y(istart:end);
            stdpr = stdpr(istart:end);
        end
        %Fit FRET fraction (fluorophore pop) v intensity
        fitmod = @(a,b,x)(a./(x+a*b*(al-1))).*(x-b).*(x>b);
        fresult = fit(x',y',fitmod,'StartPoint',[max(y),min(x)],'Weights',(1./stdpr).^2',...
            'Lower',[0,0],'Upper',[1,max(x)]);
        ab = fresult.a; bmonb = fresult.b; %a3b = 1+ab*(al-1);
        ci = confint(fresult,.68);
        %Save parameters from fits
        Nmon(ind) = bmonb;
        bpol(ind) = sum(x.*y./(ab*(1+(al-1).*y)));
        tub(ind) = sum(x);
        pfm(ind) = ab;
        FF(ind) = mean(y);
    end

%Calculate polymer amount, poly std (from fields of view) and int from each dilution series
%tuby is intensity. bpy is polymer amount.
%The values from each FOV are stored in FFr (6FOVs by 5 samples per
%dilution series).
%The values from 3229-3258 are used for finding the brightness per
%micromolar tubulin.

    if imax(loop) == 3153%matout series ending in 3153 and 3168 both belong 
        %to the same dilution series. We have to calculate the polymer
        %amount piecewise because of how the matouts were saved
        FFr = reshape(FF,3,5);
        bpolr = reshape(bpol,3,5);
        tubr = reshape(tub,3,5);
        nmonr = reshape(Nmon,3,5);
    elseif  imax(loop)==3168
        FFr(4:6,:) = reshape(FF,3,5);
        bpolr(4:6,:) = reshape(bpol,3,5);
        tubr(4:6,:) = reshape(tub,3,5);
        nmonr(4:6,:) = reshape(Nmon,3,5);
    else %here is the standard reshaping
        FFr = reshape(FF,6,5);
        bpolr = reshape(bpol,6,5);
        tubr = reshape(tub,6,5);
        nmonr = reshape(Nmon,6,5);
    end

if imax(loop)==3153 || imax(loop)==3258
    %Don't save polymer or monomer amount if we are on normalization factor
    %or first 1/2 of first dilution series series
else
    fov = 6;%Number of fields of view
    FFy(loop-1,:) = mean(FFr,1);
    bpy(loop-1,:) = mean(bpolr,1);
    tuby(loop-1,:) = mean(tubr,1);
    nmony(loop-1,:) = mean(nmonr,1);
    bperr(loop-1,:) = std(bpolr,1)/fov;
end
   %This days measurements are used to calculate normalization factors
   if imax(loop)==3258
        bpolr = reshape(bpol,6,5);
        nmonr = reshape(Nmon,6,5);
        p_norm = mean(bpolr,1);
        m_norm = mean(Nmon,1);
   end
end
%% Plots
%tuby is intensity, bpy is polymer amount
%Normalize polymer and int

%This data was taken on 6/30/16. They are on the back page of my oldest
%lab notebook
m1_day1 = [.45, .64, .63];
m2_day1 = [1.55, 1.55, 1.55];
m3_day1 = [2.88, 2.84, 2.9];
m4_day1 = [5.96, 6.07, 5.75];
m5_day1 = [9.22,9.49,9.53];

%This is data taken form 6/31/16. One high MT sample was created. Then
%Three separate m5 samples were aliquoted. The three MT series were created
%from these m5 samples
m1red = [1.16,1.15,1.19];
m1green = [1.15,1.17,1.19];
m1blue = [1.05,1.08,1.07];

m2red = [2.03,2.08,2.10];
m2green = [2.16,2.11,2.20];
m2blue = [2.04,2.05,2.04];

m3red = [3.91,3.91,3.83];
m3green = [3.73,3.74,3.73];
m3blue = [3.70,3.66,3.73];

m4red = [6.16,6.30,6.16];
m4green = [6.76,6.92,6.87];
m4blue = [6.49,6.61,6.57];

m5red = [9.92,10.20,10.28];
m5green = [10.15,10.08,9.99];
m5blue = [9.95,9.83,9.82];

nano_r = [mean(m1red),mean(m2red),mean(m3red),mean(m4red),mean(m5red)];
nano_g = [mean(m1green),mean(m2green),mean(m3green),mean(m4green),mean(m5green)];
nano_b = [mean(m1blue),mean(m2blue),mean(m3blue),mean(m4blue),mean(m5blue)];
nano_day1 = [mean(m1_day1),mean(m2_day1),mean(m3_day1),mean(m4_day1),mean(m5_day1)];

nano_day1 = nano_day1 /(80/50);
nano_r = nano_r / (180/90); 
nano_g = nano_g / (180/90); 
nano_b = nano_b / (180/90); 

for i=1:length(nano_r)
nano_err(i) = std([nano_r(i),nano_g(i),nano_b(i),nano_day1(i)],1);
end
mw = 100;
nano_drop = (1000/mw)*(mw/115)*( nano_day1 + nano_r + nano_g + nano_b ) / 4;
nano_err = (1000/mw)*(mw/115)*nano_err/sqrt(3);

x = (1000/mw)*5*[1/16, 1/8, 1/4, 1/2, 1];

norm = x(end)/(p_norm(end) + m_norm(end));
bp1 = bpy.*norm;
bp_err = std(bp1,0,1)/sqrt(3);
bp2 = mean(bp1,1);

% for i=1:3
% n_con(i,:) = bpy(i,:)/max(max(bpy));
% end

%%
figure(2);clf; hold all;
xlabel('Tubulin Concentration (\muM)');
ylabel('Measured Polymer Concentration (\muM)');
% plot(x,n_con(1,:),'+','color',[.6,.1,.1]);
% plot(x,n_con(2,:),'+','color',[.8,.1,.1]);
% plot(x,n_con(3,:),'+','color',[1,.1,.1]);
% plot(x,nano_day1,'+','color',[0.1,0.1,1]);
% plot(x,nano_r,'+','color',[0.1,0.1,.8]);
% plot(x,nano_g,'+','color',[0.1,0.1,.6]);
% plot(x,nano_b,'+','color',[0.1,0.1,.5]);
errorbar(x,bp2,bp_err,'.','color',[.8,.1,.1]);
errorbar(x,nano_drop,nano_err,'.','color',[0.1,0.1,.6]);
%plot(x,nano_drop,'o','color',[0.1,0.1,.6]);
plot(0:60,0:60,'--','color',[.6,.6,.6]);
legend('FRET-int','NanoDrop','slope=1 guideline','Location','NorthWest');




%%
% max_bpy =max(max(bpy));
% max_tuby =max(max(tuby));
% bpy = bpy ./ max_bpy;
% tuby = tuby ./ max_tuby;
% bperr = bperr./max_bpy;

%Plot polymer vs intensity
% figure(1);clf; hold all;
% axis([0 1.1 0 1.1]);
% xlabel('Tubulin Concentration');
% ylabel('Measured Polymer Amount');
% herr = errorbar(tuby,bpy,bperr,'.','Color',[.6,.6,.6]);
% yvr = reshape(bpy,15,1);%reshape(bpy,25,1);
% xvr = reshape(tuby,15,1);%reshape(tuby,25,1);
% f = fit(xvr,yvr,'poly1');
% x= 0:.01:1;
% fitvec = f.p1.*x+f.p2;
% plot(x,fitvec,'--','Color',[.6,.6,.6]);
% 
% ci95 = confint(f,.954);
% sloperr = (ci95(2,1) - ci95(1,1))/2;
% offseterr = (ci95(2,2) - ci95(1,2))/2;
% fprintf('slope is %f +/-%f \n offset is %f +/-%f',f.p1,sloperr,f.p2,offseterr);


