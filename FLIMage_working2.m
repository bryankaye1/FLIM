 clear;
set(0,'DefaulttextInterpreter','none');
addpath('Y:\Users\bkaye\cluster');

imin = 4265; imax = 4520; %split_matin = imax-imin+1;
ind =0;
tic
for i=imin:imax
    ind = ind + 1;

    try
        nstr = strcat('C:\Users\Bryan\Documents\MATLAB\data\FLIMage_cluster\',num2str(i),'.mat');
        load(nstr,'-mat','FLIMage_mat');
    catch exception
        %        try
        %             nstr = strcat('Y:\Users\bkaye\cluster\FLIMage_cluster\',num2str(i),'.mat');
        %             load(nstr,'-mat','FLIMage_mat');
        %         catch exception
        fprintf('no_FLIMage_file_detected');
        %        end
    end

    jmax = FLIMage_mat(1,1,1).jmax;
    pstart = 1+(ind-1)*jmax;
    % pend = i*(jmax/split_matin);

    for j = 1:jmax
        k = pstart+j-1;
        prest(:,k) = FLIMage_mat(1,1,j).prest;
        prestx(:,k) = FLIMage_mat(1,1,j).prestx;
        ni(k) = FLIMage_mat(1,1,j).ni;
    end
    %al = FLIMage_mat(1,1,1).w1/FLIMage_mat(1,1,1).w2;
    al = FLIMage_mat(1,1,1).w1Best/FLIMage_mat(1,1,1).w2Best;

end
toc
   

   
    %%
tic
    flimin = inf;
    for l = 1:128*128
        j = mod(l-1,128)+1;
        k = floor((l-1)/128) + 1;
        
        intensity(j,k) = ni(l); 
        [fBest,~] = transform_wf_to_f(prest(:,l), prestx(:,l),al,0,0,100, 'combine_background',10000);
        flimap(j,k)= fBest;
        if flimin>flimap(j,k)
            flimin = flimap(j,k);
        end
        
    end
toc
    whiteval = max(max(flimap));
    blackval = flimin;%min(min(flimap))+0.1;
    
    I1 = mat2gray(flimap,[blackval whiteval]);
    figure(1); clf; imshow(I1);
    
    whiteval = max(max(intensity));
    blackval = min(min(intensity));
    
    I2 = mat2gray(intensity,[blackval whiteval]);
    figure(2); clf; imshow(I2); 
    
    for l = 1:128*128
        j = mod(l-1,128)+1;
        k = floor((l-1)/128) + 1;
        x(l) =  intensity(j,k);
        y(l) = flimap(j,k);
       % stdpr(l) = flimaperr(j,k);
    end
    figure(3); clf; hold on;
    plot(x,y,'.','Markersize', 3,'Color','k');
    
    
    
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
    
    
    
    
