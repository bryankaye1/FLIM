%%change nph_counts to pixel_counts after I rerun the int_binned analysis

%close all;
clear;
Mar7_8x = [28917:28919,28921,28922];
Mar2_8x = [28386:28388];
Mar2_12x = 28391:28393;
Jan26_8x = [29472,29474,29475,29477];
Feb23_2x = [28923:28926];
alldays = [Feb23_2x,Jan26_8x,Mar2_12x,Mar2_8x,Mar7_8x];

ivec = alldays;%27252:27265;%[27252,27254:27257];
acq_time = 400;
ind = 0;
dataname_ind = {};
for i = ivec
    clear x y stdpr pixel_counts
    ind = ind + 1;
    %Load in data and images
    [~,output,flagoutput] = load_mat_data(i,1);
    %%Get FRET fraction and intensity. For now we make x(j)=ni
    ni = output(1,1,1).ni;
    al = output(1,1,1).w1Best./output(1,1,1).w2Best;
    for j = 1:length(ni)
        x(j) = ni(j)*(sum(output(1,1,j).prest'.*output(1,1,j).prestx)...
            +sum(output(1,1,j).w02est'.*output(1,1,j).w02estx));
        [y(j),stdpr(j)] = transform_wf_to_f(output(1,1,j).prest,...
            output(1,1,j).prestx,al,output(1,1,j).w02est,...
            output(1,1,j).w02estx,x(j),'dont_combine');
    end  
    intensity = x./(acq_time/length(ni));
    if i == 28917 || i == 28918
    des_ind = find(intensity>9e4);
    intensity = intensity(des_ind);
    y = y(des_ind);
    stdpr = stdpr(des_ind);
    end     
    
    [pfm(ind),nmon(ind),fresult] = fit_pf_nm(intensity,y,stdpr,al);
    ci68 = confint(fresult,.682);
    pfstd(ind) = (ci68(2,1)-ci68(1,1))/2;
    
    max_FRET(ind) = max(y);
    [~,ind_max] = max(y); 
    max_FRET_std(ind) = stdpr(ind_max);
    
    figure(ind); clf; hold on; errorbar(intensity,y,stdpr,'ko');
    xmodel = 0:max(intensity)*.001:max(intensity)*1.2;
    plot(xmodel,pf_nm_pred(al,pfm(ind),nmon(ind),xmodel),'--','Color',[0.6,0.6,0.6]);
    ti{ind} = strrep(strrep(output(1,1,1).dataname,'.sdt',''),'_',' ');
    title(ti{ind});
    
%     xintmat(ind,:) = x_int;
%     xmat(ind,:) = x;
%     ymat(ind,:) = y;
%     stdmat(ind,:) = stdpr;
%     thr_pho{ind} = find(xmat(ind,:)<10000);
    
    dataname{ind} = [output(1,1,1).pth_data,output(1,1,1).dataname];  
end

pf_est = table;
pf_est.filename = dataname';
pf_est.pf_fit = pfm';
pf_est.pf_fit_std = pfstd';
pf_est.max_FRET = max_FRET';
pf_est.max_FRET_std = max_FRET_std';
save('pf_est_table','pf_est');


% for j = 1:ind
%     figure(j); clf;
%     xmax = max(max(xintmat));
%     ymax = max(max(ymat));
%     hold on;
%     errorbar(xintmat(j,:),ymat(j,:),stdmat(j,:),'.');
%     errorbar(xintmat(j,thr_pho{j}),ymat(j,thr_pho{j}),stdmat(j,thr_pho{j}),'.','Color','red');
%     xmodel = 0:xmax*.001:xmax*1.2;
%     plot(xmodel,pf_nm_pred(al,pfm(j),Nmon(j),xmodel),'--','Color',[0.6,0.6,0.6]);
%     title(ti{j},'Interpreter','Tex');
%     xlabel('intensity (cps)');
%     ylabel('FRET fraction');
%     axis([0 xmax*1.2 0 ymax*1.2]);
% end
