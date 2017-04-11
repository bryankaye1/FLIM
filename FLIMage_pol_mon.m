clear;
%All boxcar ranges (reaches)
% low_exp = {6049:6112,6113:6176,6177:6240,6241:6304,6305:6368,6369:6432,6433:6496};
% high_exp = {6497:6560,6561:6624,6625:6688,6689:6752,6753:6816,6817:6880,6881:6944};
% low_acc = {5729:5792,5793:5856,5857:5920,5921:5984,5985:6048};
% high_acc = {5409:5472,5473:5536,5537:5600,5601:5664,5665:5728};
% superhigh_acc = {4265:4520,4777:4840,4841:4904,4905:4968,4969:5032};

%Boxcar range = 2 (reach =2, 25 pixels)
%format is matnums,colorbar fret, greybar intensity, acquistion
%number,colorbar FRET ticks

pol_cbar = [0,.5];
pol_cbar_tix = [0,.1,.2,.3,.4,.5];
pol_gbar = [0,2200];
pol_gbar_tix = 0:500:2200;
pol_low1 = {6976:7039,pol_cbar,pol_gbar,1,pol_cbar_tix,pol_gbar_tix}; %came back

sample = {29183:29310,[0,.2]};
ind =0;

for i=sample{1}
    ind = ind + 1;
    try
        nstr = strcat('/Users/bryankaye/Documents/MATLAB/data/matout/matout',num2str(i),'.mat');
        load(nstr,'-mat','FLIMage_mat');
    catch exception
        try 
            nstr = strcat('C:\Users\Bryan\Documents\MATLAB\data\FLIMage_cluster\',num2str(i),'.mat');
            load(nstr,'-mat','FLIMage_mat');
        catch
            fprintf('no_FLIMage_file_detected matnum%s\n',num2str(i));
            break;
        end
    end 
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
tic
x_grid = 0:.002:1;
count = 0;
clear fmarg flimap flimap_stat intensity

for l = 1:128*128
    j = mod(l-1,128)+1;
    k = floor((l-1)/128) + 1;
    intensity(j,k) = ni(l);
   % intensity(j,k) = (ni(l)/(sample{4}*((1+2*boxcar_range)^2)))*...
    %    (sum(prest(:,l).*prestx(:,l))+sum(w02est(:,l).*w02estx(:,l)));
    [fBest,~] = transform_wf_to_f(prest(:,l), prestx(:,l),al,...
       w02est(:,l),w02estx(:,l),600, 'n/a',10000,'mode');%,'no_histogram');

    flimap(j,k)= fBest; 
end
%count %#ok<*NOPTS>
toc%%

intensity = intensity(1+boxcar_range:end-boxcar_range,1+boxcar_range:end-boxcar_range);
flimap = flimap(1+boxcar_range:end-boxcar_range,1+boxcar_range:end-boxcar_range);

pf = .17;
pol = flimap.*intensity ./ (pf*(1+(al-1).*flimap));
mon = intensity.*(pf-flimap) ./ (pf*(1+(al-1).*flimap));
    
figure(1); clf;subplot(1,2,1);
imagesc(flimap,sample{2}); axis square;
title('FRET-fraction');
colorbar;
%set(gca,'XTickLabel','','XTick',[],'YTickLabel','','YTick',[]);
%hcb=colorbar('eastoutside');
%set(hcb,'YTick',sample{5},'YTickLabel',[]);

subplot(1,2,2); imshow(mat2gray(intensity));
title('Intensity');
colorbar;
%set(gca,'XTickLabel','','XTick',[],'YTickLabel','','YTick',[]);
%gcb = colorbar('westoutside');
%set(gcb,'Ytick',sample{6},'YTickLabel',[]);

figure(2); clf; subplot(1,2,1);
imagesc(pol); 
axis square;
title('Polymer Concentration Map');
colorbar;
%set(gca,'XTickLabel','','XTick',[],'YTickLabel','','YTick',[]);
%hcb=colorbar('eastoutside');

subplot(1,2,2); 
imagesc(mon);
axis square;
title('Monomer Concentration Map');
colorbar;
%gcb = colorbar('westoutside');


%%
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

figure(3); clf; hold on; 
plot(x,y,'.','Markersize', 3,'Color','k');
plot(x_mean,y_mean, '.', 'MarkerSize',30);    
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



