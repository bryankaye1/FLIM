%%Used to plot polymer FF vs pixel group for spindle area FRET fraction
%%measurements.

function [ymat,varmat,imat,polmat,monmat,mask_distance] = FF_vs_dist_int(ivec,...
    average_spindles,show_spindle,pf,showmon)

ind = 0;
load('/Users/bryankaye/Documents/MATLAB/FLIM/pf_est_table');
for i = ivec
    clear x y distance stdpr
    ind = ind+1;
    [~,output] = load_mat_data(i,'local',1);
    load(['/Users/bryankaye/Documents/MATLAB/data/matin/matin',...
        num2str(i),'.mat']);
    if seg_results.input_params.scan_mag == 8
        maxlen = 17;%20;
        minlen = -5;%-5.75;
    elseif seg_results.input_params.scan_mag == 12.8
        maxlen = 8;%10;
        minlen = -5;%-5.75;
    elseif seg_results.input_params.scan_mag == 2
        maxlen = 100;
        minlen = -5;
    end
    
    ind_maxlen_range = find(maxlen > seg_results.mask_distance);
    if average_spindles
        ind_minlen_range = find(seg_results.mask_distance > minlen);
        des_ind = intersect(ind_maxlen_range,ind_minlen_range);
    else
        des_ind = ind_maxlen_range;
    end
    
    intensity = output(1,1,1).ni;
    al = output(1,1,1).w1Best./output(1,1,1).w2Best;
    for j = 1:length(intensity)
        %num_pho(j) = intensity(j)*sum(sum(seg_results.int_masks(:,:,j)));
        num_pho = sum(output(1,1,j).datahis);
        [y(j),stdpr(j)] = transform_wf_to_f(output(1,1,j).prest,...
            output(1,1,j).prestx,al,output(1,1,j).w02est,...
            output(1,1,j).w02estx,501,'dont_combine'); %#ok<*SAGROW>
    end
    
    mask_distance = seg_results.mask_distance(des_ind);
    intensity = intensity(des_ind);
    y = y(des_ind);
    stdpr = stdpr(des_ind);
    [pfm(ind),nmon(ind)] = fit_pf_nm(intensity,y,stdpr,al);
    
    figure(ind); clf; hold on; errorbar(intensity,y,stdpr,'ko');
    xmodel = 0:max(intensity)*.001:max(intensity)*1.2;
    plot(xmodel,pf_nm_pred(al,pfm(ind),nmon(ind),xmodel),'--','Color',[0.6,0.6,0.6]);
    xlabel('intensity (au)');
    ylabel('FRET fraction');
    
    ti{ind} = [output(1,1,1).pth_data(end-10:end),...
        strrep(strrep(output(1,1,1).dataname,'.sdt',''),'_',' ')];
    title(ti{ind});
    
    pol = y.*intensity ./ (pf*(1+(al-1).*y));
    mon = intensity.*(pf-y) ./ (pf*(1+(al-1).*y));
    mon = mon/max(mon);
    pol = pol/ max(pol);
    
    des_row = find(ismember(pf_est.filename,...
        [output(1,1,1).pth_data,output(1,1,1).dataname]));
    % pf = max([pf_est.pf_fit(des_row)+pf_est.pf_fit_std(des_row),...
    %      pf_est.max_FRET(des_row)+pf_est.max_FRET_std(des_row)]);
    
    pf = max([pf_est.pf_fit(des_row),pf_est.max_FRET(des_row)]);
    %pf = pf_est.max_FRET(des_row);
    
    % pf = pf_est.max_FRET(des_row)+pf_est.max_FRET_std(des_row);
    
    
    
    if average_spindles
        imat(ind,:) = intensity;
        ymat(ind,:) = y;
        polmat(ind,:) = pol;
        monmat(ind,:) = mon;
        varmat(ind,:) = stdpr.*stdpr;
    end
    
end


if average_spindles
    y_ave_weighted = sum(ymat./varmat,1);
    norm = sum(1./varmat,1);
    y_ave = y_ave_weighted./norm;
    %y_ave = mean(ymat,1);
    i_ave = mean(imat,1);
    %pol_ave = mean(polmat,1);
    pol_ave = y_ave.*i_ave ./ (pf*(1+(al-1).*y_ave));
    pol_ave = pol_ave/max(pol_ave);
    mon_ave = mean(monmat,1);
    std_ave = sqrt(mean(varmat,1)/(ind));
    
    ti = [];
    if showmon
        figure(ind+1); clf;
        plot(mask_distance,mon_ave,'bo');
        title('monomer average');
        xlabel('distance (microns)');
        ylabel('Monomer Concentration (au)');
        axis([min(mask_distance) maxlen 0 1]);
    else
        FFdist_plots(ind+1,mask_distance,maxlen,pol_ave,i_ave,y_ave,...
            std_ave,ti,0,seg_results);
    end
end
end




%%




%     des_row = find(ismember(pf_est.filename,...
%         [output(1,1,1).pth_data,output(1,1,1).dataname]));
%     pf = max([pf_est.pf_fit(des_row)+pf_est.pf_fit_std(des_row),...
%         pf_est.max_FRET(des_row)+pf_est.max_FRET_std(des_row)]);
%
%      pf = max([pf_est.pf_fit(des_row),pf_est.max_FRET(des_row)]);
%     pf = pf_est.max_FRET(des_row)+pf_est.max_FRET_std(des_row);
%
%     if isempty(pf)
%         pf = 0.15;
%     end



%    num_pixels = squeeze(sum(sum(seg_results.int_masks,1),2))';

% %     ind_flat = find(10 > seg_results.mask_distance & seg_results.mask_distance > 5);
% %     int_flat = mean(intensity(ind_flat));
% %     fret_flat = mean(y(ind_flat));
%     fret_peak =max(y);
%     int_peak = max(intensity);
% %     scale =  (fret_flat - fret_peak) / (int_flat - int_peak);
% %     offset = ( fret_peak * int_flat - fret_flat * int_peak ) / (int_flat - int_peak);
% %     intensity_norm = intensity *scale + offset;
%     intensity_norm = intensity * fret_peak/int_peak;
%     plot(mask_distance,intensity_norm,'ro');
%     legend('FRET-fraction','Intensity');
%   figure(ind+10); imshow(mat2gray(seg_results.image_stack));


%clear int_image
% int_image{1} = seg_results.first_image;
% int_image{2} = seg_results.middle_image;
% int_image{3} = seg_results.last_image;
% image_stack = seg_results.image_stack;
%
% masks = seg_results.int_masks;
%
%
% h2 = figure(2); clf;
% subplot(2,2,1); imshow(mat2gray(int_image{1}),'InitialMagnification','fit');
% title('first intensity image');
% subplot(2,2,2); imshow(mat2gray(int_image{2}),'InitialMagnification','fit');
% title('middle intensity image');
% subplot(2,2,3); imshow(mat2gray(int_image{end}),'InitialMagnification','fit');
% title('last intensity image');
% subplot(2,2,4); imshow(mat2gray(image_stack),'InitialMagnification','fit');
% title('registered intensity image');
%
%
% h3 = figure(3); clf;
% B= imoverlay(mat2gray(image_stack),max(masks(:,:,8),masks(:,:,end)));
% imshow(B,'InitialMagnification','fit');
% %imagesc(int_final)
% %axis equal
% %hold on;
% title('aligned registered image and first and last mask');
