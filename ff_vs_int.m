%%change nph_counts to pixel_counts after I rerun the int_binned analysis

%close all;
clear;
ivec = 28377:28381;%27252:27265;%[27252,27254:27257];
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
    %al = .45;
    for j = 1:length(ni)
        x(j) = ni(j)*(sum(output(1,1,j).prest'.*output(1,1,j).prestx)...
            +sum(output(1,1,j).w02est'.*output(1,1,j).w02estx));
        [y(j),stdpr(j)] = transform_wf_to_f(output(1,1,j).prest,...
            output(1,1,j).prestx,al,output(1,1,j).w02est,...
            output(1,1,j).w02estx,x(j),'dont_combine');
    end
    
    
    %num_pixels = squeeze(sum(sum(int_masks,1),2))';
    intensity = x./(acq_time/length(ni));   
    
    [pfm(ind),nmon(ind),fresult] = fit_pf_nm(intensity,y,stdpr,al);
    ci68 = confint(fresult,.682);
    pfstd(ind) = (ci68(2,1)-ci68(1,1))/2;
    
    figure(ind); clf; hold on; errorbar(intensity,y,stdpr,'ko');
    xmodel = 0:max(intensity)*.001:max(intensity)*1.2;
    plot(xmodel,pf_nm_pred(al,pfm(ind),nmon(ind),xmodel),'--','Color',[0.6,0.6,0.6]);
    ti{ind} = strrep(output(1,1,1).dataname,'.sdt',''); title(ti{ind});
    
%     xintmat(ind,:) = x_int;
%     xmat(ind,:) = x;
%     ymat(ind,:) = y;
%     stdmat(ind,:) = stdpr;
%     thr_pho{ind} = find(xmat(ind,:)<10000);
    
end

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
