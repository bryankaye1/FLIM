%%Used to plot polymer series plots

clear;
figure(1);clf;
imin = [3139,3154,3169,3199];%,3229];%
imax = [3153,3168,3198,3228];%,3258]; %
for loop = 1:4

    ind = 0; %index reset to zero for each dilution series
    %%Load in data for one dilution series, fit FF v int,
    %%and calculate polymer amount given fit parameters
    for i = imin(loop):imax(loop)  %Load in matout files
        clear output intb ffP cintpr
        [~,output] = load_mat_data(i);
        ind =ind+1;
        dataname = output.dataname;
        %%Set fit paramters for FF v int over image
        ni = output(1,1,1).ni; 
        al = output(1,1,1).w1Best./output(1,1,1).w2Best;
        ep = output(1,1,1).sinti; 
        %%Get FRET fraction (fluoro population) and intensity
        for j = 1:length(output)
            x(j) = (ni(j)/ep(j))*(sum(output(1,1,j).prest'.*output(1,1,j).prestx)...
               +sum(output(1,1,j).w02est'.*output(1,1,j).w02estx));
            [y(j),stdpr(j)] = transform_wf_to_f(output(1,1,j).prest,...
                output(1,1,j).prestx,al,output(1,1,j).w02est,output(1,1,j).w02estx,x(j));
        end
        %Remove pixel groups with zero photons
        if min(x)==0
            istart = find(x==0,1,'last')+1;
            x = x(istart:end);
            y = y(istart:end);
            stdpr = stdpr(istart:end);
        end
        %Fit FRET fraction (fluoro pop) v intensity
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
        pfm_name{ind} = {output(1,1,1).dataname,ab};
    end

    %%%%Calculate polymer amount, poly std (from fields of view) and 
    %%%%int from each dilution series
        fov = 6;%Number of fields of view
    %tuby is intensity
    %bpy is polymer amount

    if imax(loop) == 3153%matout series ending in 3153 and 3168 both belong 
        %to the same dilution series. We have to calculate the polymer
        %amoun piecewise because of how the matouts were saved
        bpolr = reshape(bpol,3,5);
        tubr = reshape(tub,3,5);
    elseif  imax(loop)==3168
        bpolr(4:6,:) = reshape(bpol,3,5);
        tubr(4:6,:) = reshape(tub,3,5);
        bpy(loop-1,:) = mean(bpolr,1);
        tuby(loop-1,:) = mean(tubr,1);
        bperr(loop-1,:)=std(bpolr,1)/fov; 
        
    elseif imax(loop)==3198
        bpolr = reshape(bpol,6,5);
        tubr = reshape(tub,6,5);
        bpy(loop-1,:) = mean(bpolr,1);
        tuby(loop-1,:) = mean(tubr,1);
        bperr(loop-1,:)=std(bpolr,1)/fov;
    elseif imax(loop)==3228
        bpolr = reshape(bpol,6,5);
        tubr = reshape(tub,6,5);
        bpy(loop-1,:) = mean(bpolr,1);
        tuby(loop-1,:) = mean(tubr,1);
        bperr(loop-1,:)=std(bpolr,1)/fov;
    end
junk_var = 0;

end
%% Plots
%tuby is intensity, bpy is polymer amount
%Normalize polymer and int
max_bpy =max(max(bpy));
max_tuby =max(max(tuby));
bpy = bpy ./ max_bpy;
tuby = tuby ./ max_tuby;
bperr = bperr./max_bpy;

%Plot polymer vs intensity
figure(1);clf; hold all;
axis([0 1.1 0 1.1]);
xlabel('Tubulin Concentration');
ylabel('Measured Polymer Amount');
herr = errorbar(tuby,bpy,bperr,'.','Color',[.6,.6,.6]);
yvr = reshape(bpy,15,1);%reshape(bpy,25,1);
xvr = reshape(tuby,15,1);%reshape(tuby,25,1);
f = fit(xvr,yvr,'poly1');
x= 0:.01:1;
fitvec = f.p1.*x+f.p2;
plot(x,fitvec,'--','Color',[.6,.6,.6]);

ci95 = confint(f,.954);
sloperr = (ci95(2,1) - ci95(1,1))/2;
offseterr = (ci95(2,2) - ci95(1,2))/2;
fprintf('slope is %f +/-%f \n offset is %f +/-%f',f.p1,sloperr,f.p2,offseterr);


