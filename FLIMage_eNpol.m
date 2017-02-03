clear;
plotFF = 0;
%All boxcar ranges (reaches)
% low_exp = {6049:6112,6113:6176,6177:6240,6241:6304,6305:6368,6369:6432,6433:6496};
% high_exp = {6497:6560,6561:6624,6625:6688,6689:6752,6753:6816,6817:6880,6881:6944};
% low_acc = {5729:5792,5793:5856,5857:5920,5921:5984,5985:6048};
% high_acc = {5409:5472,5473:5536,5537:5600,5601:5664,5665:5728};
% superhigh_acc = {4265:4520,4777:4840,4841:4904,4905:4968,4969:5032};

%Boxcar range = 2 (reach =2, 25 pixels)
%format is matnums,colorbar fret, greybar intensity, acquistion
%number,colorbar FRET ticks
low_exp = {6177:6240,[0,.4],[0,660],1,[0,.1,.2,.3,.4],0:100:600};
high_exp = {6625:6688,[0,.4],[0,660],1,[0,.1,.2,.3,.4],0:100:600};
low_acc = {5857:5920,[0,.15],[0,300],10,[0,.05,.1,.15,.2],[0,50,100,150,200,250]};
high_acc = {5537:5600,[0,.15],[0,300],10,[0,.05,.1,.15,.2],[0,50,100,150,200,250]};
superhigh_acc = {4841:4904,[0,.2],[0,277],10,[0,.05,.1,.15,.2],[0,50,100,150,200,250]}; %index
notax_superhigh_acc = {26346:26473,[0,.2],[0,277],10,[0,.05,.1,.15,.2],[0,50,100,150,200,250]};
%25956:26083 for 2p5A, which had a higher pf when taxol was added

no_taxol_too_high_superhigh_acc_ = {25796:25923,[0,.2],[0,277],10,[0,.05,.1,.15,.2],[0,50,100,150,200,250]};

feb27_reach0 = {27266:27393,0,0,1};
m2 = {27970:28033,0,0,1};
m4 = {27842:27905,0,0,1};
m5 = {27906:27969,0,0,1};

m3r0 = {27778:27841,0,0,1};
m3 = {27650:27777,0,0,1};
m3r2 = {27522:27649,0,0,1};



pol_cbar = [0,.5];
pol_cbar_tix = [0,.1,.2,.3,.4,.5];
pol_gbar = [0,2200];
pol_gbar_tix = 0:500:2200;

pol_low1 = {6976:7039,pol_cbar,pol_gbar,1,pol_cbar_tix,pol_gbar_tix}; %came back
pol_low2 = {7040:7103,pol_cbar,pol_gbar,1,pol_cbar_tix,pol_gbar_tix};
pol_low3 = {7104:7167,pol_cbar,pol_gbar,1,pol_cbar_tix,pol_gbar_tix};

pol_hi1 = {7168:7231,pol_cbar,pol_gbar,1,pol_cbar_tix,pol_gbar_tix};
pol_hi2 = {7232:7295,pol_cbar,pol_gbar,1,pol_cbar_tix,pol_gbar_tix};
pol_hi3 = {7296:7359,pol_cbar,pol_gbar,1,pol_cbar_tix,pol_gbar_tix};
pol_hi4 = {7360:7423,pol_cbar,pol_gbar,1,pol_cbar_tix,pol_gbar_tix};
pol_hi5 = {7424:7487,pol_cbar,pol_gbar,1,pol_cbar_tix,pol_gbar_tix};
pol_hi6 = {7488:7551,pol_cbar,pol_gbar,1,pol_cbar_tix,pol_gbar_tix};
pol_cell = {pol_low1,pol_low2,pol_low3,pol_hi1,pol_hi2,pol_hi3,pol_hi4,...
    pol_hi5,pol_hi6};

spindle_jan5 = {27090:27217};

plotind = 0;
%for j_outer = 1:length(pol_cell)

sample = m4;%pol_cell{j_outer};


%for keyind = superhigh_acc
%    clear prest prestx
ind =0;
%for i=keyind{1}
tic
for i=sample{1}
    ind = ind + 1;
    
    nstr = strcat('/Users/bryankaye/Documents/MATLAB/data/matout/matout',num2str(i),'.mat');
    load(nstr,'-mat');
    
    jmax = FLIMage_mat(1,1,1).jmax;
    pstart = 1+(ind-1)*jmax;
    % pend = i*(jmax/split_matin);
    
    for j = 1:jmax
        k = pstart+j-1;
        prest(:,k) = FLIMage_mat(1,1,j).prest;
        prestx(:,k) = FLIMage_mat(1,1,j).prestx;
        w02est(:,k) = FLIMage_mat(1,1,j).w02est;
        w02estx(:,k) = FLIMage_mat(1,1,j).w02estx;
        % ni(k) = FLIMage_mat(1,j,1).ni;
    end
    %al = FLIMage_mat(1,1,1).w1/FLIMage_mat(1,1,1).w2;
    al = FLIMage_mat(1,1,1).w1Best/FLIMage_mat(1,1,1).w2Best;
    
end
ni = FLIMage_mat(1,1,1).ni;
boxcar_range = FLIMage_mat(1,1,1).reach;

%%

%flimin = inf;
x_grid = 0:.002:1;
count = 0;
clear fmarg flimap flimap_stat intensity

for l = 1:128*128
    j = mod(l-1,128)+1;
    k = floor((l-1)/128) + 1;
    
    intensity(j,k) = (ni(l)/(sample{4}*((1+2*boxcar_range)^2)))*...
        (sum(prest(:,l).*prestx(:,l))+sum(w02est(:,l).*w02estx(:,l)));
    [fBest,~] = transform_wf_to_f(prest(:,l), prestx(:,l),al,...
        w02est(:,l),w02estx(:,l),100, 'n/a',2000,'no_histogram');
    
    flimap(j,k)= fBest;
    
end

intensity = intensity(1+boxcar_range:end-boxcar_range,1+boxcar_range:end-boxcar_range);
flimap = flimap(1+boxcar_range:end-boxcar_range,1+boxcar_range:end-boxcar_range);

 %%

pf = 0.49;
norm_FOV = 2.680e-6;
%divided by 10 for conversation from uM to mg/ml. 
%divide by total number of pixels to convert from FOV to per pixel norm
norm = norm_FOV * (126*126)/10; 

polmap = norm*intensity.*flimap ./ (pf*(1+(al-1).*flimap));
figure; clf; imagesc(polmap); axis square;
%set(gca,'XTickLabel','','XTick',[],'YTickLabel','','YTick',[]);
gcb = colorbar('westoutside');
%set(gcb,'Ytick',sample{6},'YTickLabel',[]); drawnow;

toc




%%







if plotFF
    for l = 1:(128-2*boxcar_range)*(128-2*boxcar_range)
        j = mod(l-1,(128-2*boxcar_range))+1;
        k = floor((l-1)/(128-2*boxcar_range)) + 1;
        x(l) =  intensity(j,k);
        y(l) = flimap(j,k);
        % stdpr(l) = flimaperr(j,k);
    end
    num_int_bins = 16;
    intbin = min(x)+(0:num_int_bins)*(max(x)-min(x))/num_int_bins;
    x_mean = [];
    y_mean = [];
    yp_mean =[];
    yp_mode =[];
    for k = 1:num_int_bins
        x_total = 0;
        y_total = 0;
        yp_total = ones(length(x_grid),1);
        count =0;
        for l = 1:length(x)
            if intbin(k) < x(l) && intbin(k+1)> x(l)
                count = count  +1;
                x_total = x(l)+x_total;
                y_total = y(l)+y_total;
                
                %yp_total = f_marg(:,l).*yp_total;
                if sum(isinf(yp_total)) || sum(isnan(yp_total))
                    l
                end
                %[fBest,~,fest,festx] = transform_wf_to_f_suppress_isnan(tempf.output(1,1,l).prest,...
                %          tempf.output(1,1,l).prestx,al,tempf.output(1,1,l).w02est,...
                %           tempf.output(1,1,l).w02estx,100);
                %yp_total= fest .* yp_total;
            end
        end
        if count > 0
            x_mean =[x_mean x_total/count];
            y_mean = [y_mean y_total/count];
            
            if sum(yp_total)==0
                yp_norm = 0;
                ind2 = 1;
            else
                yp_norm = yp_total/sum(yp_total);
                [~,ind2] = max(yp_total);
            end
            
            yp_mean = [yp_mean sum(yp_norm.*x_grid')];
            yp_mode = [yp_mode x_grid(ind2)];
            %yp_mean = [yp_mean ypmax];
        end
    end
    figure(plotind*3+3); clf; hold on;
    plot(x,y,'.','Markersize', 3,'Color','k');
    
    if sample{1}(1)==4841
        set(gca,'XTick',0:100:300,'Ytick',0:0.04:.2,'YTickLabel',[],'XTickLabel',[]); %was 'Ytick',0:0.03:.12
        axis([15 350 0 .2]); %Was axis([15 350 0 .12]); in original submission
    else
        plot(x_mean,y_mean, '.', 'MarkerSize',30);
    end
    plotind = plotind+1;
end




%end
%ti = sprintf('Fret Fraction vs Intensity: boxcar range is %s',num2str(boxcar_range)); title(ti);
%xlabel('Number of photons');
%ylabel('FRET population-fraction'); drawnow;





%%

%%
%     clear flimap
%     al = tempf.output(1,1,1).w1Best/tempf.output(1,1,1).w2Best;
%     flimin = pi;
%     tic
%     for l = 1:128*128
%         j = mod(l-1,128)+1;
%         k = floor((l-1)/128) + 1;
%
%         intensity(j,k) = sum(tempf.output(1,1,l).datahis);
%         %         [fBest,stdpr] = transform_wf_to_f_suppress_isnan(tempf.output(1,1,l).prest,...
%         %             tempf.output(1,1,l).prestx,al,tempf.output(1,1,l).w02est,...
%         %             tempf.output(1,1,l).w02estx,100);
%         %         if (fBest-5*stdpr)<0.02
%         %             flimap(j,k) = 0;
%
%         [~,flimaperr(j,k),~] = findci(tempf.output(1,1,l).prest, tempf.output(1,1,l).prestx, 0.68,'error_size');
%         [~,rel_heights,~] = findci(tempf.output(1,1,l).prest, tempf.output(1,1,l).prestx, 0.68,'dist_endpoints');
%         if 1==0%rel_heights < 5
%             flimap(j,k)=0;
%
%             %          if dist_zero <0.05
%             %              flimap(j,k)=0;
%
%             % %          if intensity(j,k) < 1000
%             % %              flimap(j,k) = 0;
%         else
%                          [fBest,~] = transform_wf_to_f(tempf.output(1,1,l).prest,...
%                              tempf.output(1,1,l).prestx,al,tempf.output(1,1,l).w02est,...
%                              tempf.output(1,1,l).w02estx,100, 'combine_background');
%             %flimap(j,k)= tempf.output(1,1,l).prBest;%tempf.output(1,1,l).prBest/(tempf.output(1,1,l).prBest+al*tempf.output(1,1,l).w02Best);
%             flimap(j,k)= fBest;
%             if flimin>flimap(j,k)
%                 flimin = flimap(j,k);
%             end
%         end
%     end
% toc
%     whiteval = max(max(flimap));
%     blackval = flimin;%min(min(flimap))+0.1;
%
%     I1 = mat2gray(flimap,[blackval whiteval]);
%     figure(1); clf; imshow(I1); title(tempf.output(1,1,1).dataname);
%
%     whiteval = max(max(intensity));
%     blackval = min(min(intensity));
%
%     I2 = mat2gray(intensity,[blackval whiteval]);
%     figure(2); clf; imshow(I2); title(tempf.output(1,1,1).dataname);
%
%
%     %%
%     size_flimap = size(flimap);
%     clear x y stdpr
%
%     for l = 1:128*128
%         j = mod(l-1,128)+1;
%         k = floor((l-1)/128) + 1;
%         x(l) =  intensity(j,k);
%         y(l) = flimap(j,k);
%         stdpr(l) = flimaperr(j,k);
%     end
%     figure(3); clf; hold on;
%     plot(x,y,'.','Markersize', 3,'Color','k');
% %%
%     num_int_bins = 16;
%     intbin = min(x)+(0:num_int_bins)*(max(x)-min(x))/num_int_bins;
%     x_mean = [];
%     y_mean = [];
%     yp_mean =[];
%     for k = 1:num_int_bins
%         x_total = 0;
%         y_total = 0;
%         yp_total = 1;
%         count =0;
%         for l = 1:length(x)
%             if intbin(k) < x(l) && intbin(k+1)> x(l)
%                 count = count  +1;
%                 x_total = x(l)+x_total;
%                 y_total = y(l)+y_total;
%                 %[fBest,~,fest,festx] = transform_wf_to_f_suppress_isnan(tempf.output(1,1,l).prest,...
%                 %          tempf.output(1,1,l).prestx,al,tempf.output(1,1,l).w02est,...
%                 %           tempf.output(1,1,l).w02estx,100);
%                 %yp_total= fest .* yp_total;
%             end
%         end
%         if count > 0
%             x_mean =[x_mean x_total/count];
%             y_mean = [y_mean y_total/count];
%             %[~,ind] = max(yp_total);
%             %ypmax = festx(ind);
%             %yp_mean = [yp_mean ypmax];
%         end
%     end
%
%        figure(3); hold on;
%       plot(x_mean,y_mean, 'o');
%
%
%
%


%h2=errorbar(x,y,stdpr,'.','Color',[.6,.6,.6]);


%     for j = 1:length(output)
%         x(j) = (ni(j)/ep(j))*(sum(output(1,1,j).prest'.*output(1,1,j).prestx)...
%             +sum(output(1,1,j).w02est'.*output(1,1,j).w02estx));
%         [y(j),stdpr(j)] = transform_wf_to_f(output(1,1,j).prest,...
%             output(1,1,j).prestx,al,output(1,1,j).w02est,output(1,1,j).w02estx,x(j));
%     end
%Remove pixel groups with zero photons
%     if min(x)==0
%         istart = find(x==0,1,'last')+1;
%         x = x(istart:end);
%         y = y(istart:end);
%         cintpr = cintpr(istart:end);
%     end
%     %Fit FRET fraction (fluoro pop) v intensity
%     fitmod = @(a,b,x)(a./(x+a*b*(al-1))).*(x-b).*(x>b);
%     fresult = fit((x)',y',fitmod,'StartPoint',[max(y),min(x)],'Weights',(1./stdpr).^2',...
%         'Lower',[0,0],'Upper',[1,max(x)]);
%     ab = fresult.a; bmonb = fresult.b;
%     ci95 = confint(fresult,.954);
%     ci68 = confint(fresult,.682);
%     xmodel = 1000:1000:8*10^6;
%     modelb = fitmod(ab,bmonb,xmodel);



