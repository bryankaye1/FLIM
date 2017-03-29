%%Used to plot polymer FF vs pixel group for spindle area FRET fraction
%%measurements.

clear;
ivec = 28519:28524;%28330:28332; %
ind = 0;
average_spindles = 0;


for i = ivec
    clear x y distance stdpr
    ind = ind+1;
    [~,output] = load_mat_data(i,'local',1);
    load(['/Users/bryankaye/Documents/MATLAB/data/matin/matin',...
        num2str(i),'.mat']);
    int_masks = seg_results.int_masks;
    mask_distance = seg_results.mask_distance;
    ni = output(1,1,1).ni;
    al = output(1,1,1).w1Best./output(1,1,1).w2Best;

  %  figure(ind+10); imshow(mat2gray(seg_results.image_stack));
    %%Get FRET fraction (fluoro population) and intensity
    
    for j = 1:length(output)
        x(j) = ni(j)*(sum(output(1,1,j).prest'.*output(1,1,j).prestx)...
            +sum(output(1,1,j).w02est'.*output(1,1,j).w02estx)); %#ok<SAGROW>
        [y(j),stdpr(j)] = transform_wf_to_f(output(1,1,j).prest,...
            output(1,1,j).prestx,al,output(1,1,j).w02est,...
            output(1,1,j).w02estx,x(j),'dont_combine'); %#ok<*SAGROW>
    end 
    
    if seg_results.input_params.scan_mag==8
    maxlen = 20;
    elseif seg_results.input_params.scan_mag==12.8
        maxlen = 10;
    elseif seg_results.input_params.scan_mag== 2
        maxlen = 100;
    end
    pf = .15;
    

    
    thr_pho = find(ni > 0);%0.5e5/20);
    num_pixels = squeeze(sum(sum(int_masks,1),2))';
    intensity = ni;%./num_pixels; %try ni here  

    ind_maxlen_range = find(maxlen > seg_results.mask_distance);
    intensity = intensity(ind_maxlen_range);
    y = y(ind_maxlen_range);
    stdpr = stdpr(ind_maxlen_range);
    mask_distance = mask_distance(ind_maxlen_range);
    
    pol = y.*intensity ./ (pf*(1+(al-1).*y));
    mon = intensity.*(pf-y) ./ (pf*(1+(al-1).*y));
    mon = mon/max(mon);
    pol = pol/ max(pol);
    
    figure(ind);clf; 
    
    subplot(1,3,2);
    plot(mask_distance,mon,'go');
    xlabel('distance (microns)');
    ylabel('Monomer Concentration (au)');
    axis([min(mask_distance) maxlen 0 1]);
    
    subplot(1,3,3);
    plot(mask_distance,pol,'go');
    xlabel('distance (microns)');
    ylabel('Polymer Concentration (au)');
    axis([min(mask_distance) maxlen 0 1]);
    
    
    subplot(1,3,1);
    % plot(mask_distance(thr_pho),y(thr_pho),'bo')
    errorbar(mask_distance,y,stdpr,'bo');
    %errorbar(distance,y,stdpr);
    xlabel('distance (microns)');
    ylabel('FRET fraction');
    title(strrep(output(1,1,1).dataname,'_',' '));
    axis([min(mask_distance) maxlen 0 .18]);
    % axis([-2 110 0 .13]);
    
    hold on;
    
    
    ind_flat = find(10 > seg_results.mask_distance & seg_results.mask_distance > 5);
    int_flat = mean(intensity(ind_flat));
    fret_flat = mean(y(ind_flat));
    fret_peak =max(y);
    int_peak = max(intensity); 
    scale =  (fret_flat - fret_peak) / (int_flat - int_peak);
    offset = ( fret_peak * int_flat - fret_flat * int_peak ) / (int_flat - int_peak);     

    intensity_norm = intensity *scale + offset;
    
    plot(mask_distance,intensity_norm,'ro');
    legend('FRET-fraction','Intensity');
    

    
    
%     subplot(1,2,2); %errorbar(distance,intensity,sqrt(intensity));
%     plot(mask_distance(thr_pho),intensity(thr_pho),'ro');
%     xlabel('distance');
%     ylabel('intensity');
%     title(strrep(output(1,1,1).dataname,'_',' '));
    
    if average_spindles
        nmat{ind} = ni;
        xmat{ind} = mask_distance;
        ymat{ind} = y; 
        len = length(x);
        if len > maxlen
        maxlen = len; 
        max_len_ind = ind;
        end
    end
end

if average_spindles
    x_sq = nan(ind,max(maxlen));
    y_sq = nan(ind,max(maxlen));
    n_sq = zeros(ind,max(maxlen));
    
    for k = 1:ind
        x_sq(k,1:length(xmat{k})) = xmat{k};
        y_sq(k,1:length(ymat{k})) = ymat{k};
        n_sq(k,1:length(ymat{k})) = nmat{k};
    end
    
    x_ave = nanmean(x_sq,1);
    y_ave = nanmean(y_sq,1);
    n_ave = sum(n_sq,1);
    thr_phoa = find(n_ave > 1e6);
    
    figure(ind+1); plot(x_ave(thr_phoa),y_ave(thr_phoa),'bo');
    xlabel('distance from spindle edge');
    ylabel('FRET fraction');
end
%%
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
