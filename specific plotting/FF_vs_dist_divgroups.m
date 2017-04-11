%%Used to plot polymer FF vs pixel group for spindle area FRET fraction
%%measurements.

clear;
zoom8x_Mar7 = [28519:28521,28523,28524];
zoon8x_Mar2 = 28519:28512;
zoom8x_Jan26 = 28491:28493;
zoom12x = 28515:28517;
zoom8x = [zoom8x_Mar7,zoon8x_Mar2,zoom8x_Jan26];
zoom2x = 28506:28509;

zoom12x_donor = 28518;
zoom8x_donor = [28522,28513:28514,28487,28527];%28527 is from 2-28, where there was no positive control
zoom2x_donor = [28525,28526];

quad_div = 29459;

ivec = quad_div;%28330:28332; %
ind = 0;
average_spindles = 1;
show_spindle = 1;


for i = ivec
    clear x y distance stdpr
    ind = ind+1;
    [~,output] = load_mat_data(i,'local',1);
    load(['/Users/bryankaye/Documents/MATLAB/data/matin/matin',...
        num2str(i),'.mat']);
    if seg_results.input_params.scan_mag == 8
        maxlen = 20;
        minlen = -5.75;
    elseif seg_results.input_params.scan_mag == 12.8
        maxlen = 10;
        minlen = -5.75;
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
        [y(j),stdpr(j)] = transform_wf_to_f(output(1,1,j).prest,...
            output(1,1,j).prestx,al,output(1,1,j).w02est,...
            output(1,1,j).w02estx,501,'dont_combine'); %#ok<*SAGROW>
    end 
    
                pf = .24;
    switch isfield(seg_results,'div_groups')
        case 0
            mask_distance = seg_results.mask_distance(des_ind);
            intensity = intensity(des_ind);
            y = y(des_ind);
            stdpr = stdpr(des_ind);
            
            pol = y.*intensity ./ (pf*(1+(al-1).*y));
            mon = intensity.*(pf-y) ./ (pf*(1+(al-1).*y));
            mon = mon/max(mon);
            pol = pol/ max(pol);
            
            ti = strrep(output(1,1,1).dataname,'_',' ');
            %fprintf('%s is %3.0f\n',ti,i);
            FFdist_plots(ind,mask_distance,maxlen,mon,pol,intensity,y,stdpr,ti,...
                show_spindle,seg_results);
            
            num_pho = intensity.*sum(sum(seg_results.int_masks));
            if sum(num_pho<500)
                fprintf('Photons count less than 500!'); pause;
            end
            
            if average_spindles
                imat(ind,:) = intensity;
                ymat(ind,:) = y;
                polmat(ind,:) = pol;
                monmat(ind,:) = mon;
                varmat(ind,:) = stdpr.*stdpr;
            end
            
        case 1
            mask_distance = seg_results.mask_distance(des_ind);
            for dg=1:length(seg_results.div_groups)
                dg_ind = des_ind+(dg-1)*length(seg_results.div_groups{1});                
                intensity_dg = intensity(dg_ind);
                y_dg = y(dg_ind);
                stdpr_dg = stdpr(dg_ind);
                
                pol = y_dg.*intensity_dg ./ (pf*(1+(al-1).*y_dg));
                mon = intensity_dg.*(pf-y_dg) ./ (pf*(1+(al-1).*y_dg));
                mon = mon/max(mon);
                pol = pol/ max(pol);
                
                ti = strrep(output(1,1,1).dataname,'_',' ');
                %fprintf('%s is %3.0f\n',ti,i);
                FFdist_plots(ind-1+dg,mask_distance,maxlen,mon,pol,...
                    intensity_dg,y_dg,stdpr_dg,ti,show_spindle,seg_results);
                
                figure(ind-1+dg); subplot(1,4,4);
                B= imoverlay(mat2gray(seg_results.image_stack),+...
                seg_results.div_masks{dg}(:,:,15)+seg_results.div_masks{dg}(:,:,20));
                imshow(B,'InitialMagnification','fit');
                num_pho = intensity(1:length(seg_results.int_masks))...
                    .*sum(sum(seg_results.int_masks));
                if sum(num_pho<500)
                    fprintf('Photons count less than 500!'); pause;
                end
                
                if average_spindles
                    imat(dg,:) = intensity_dg;
                    ymat(dg,:) = y_dg;
                    polmat(dg,:) = pol;
                    monmat(dg,:) = mon;
                    varmat(dg,:) = stdpr_dg.*stdpr_dg;
                end  
            end            
    end   
end


if average_spindles    
    
    y_ave = mean(ymat,1);
    i_ave = mean(imat,1);
    pol_ave = mean(polmat,1);
    mon_ave = mean(monmat,1);
    std_ave = sqrt(mean(varmat,1)/(ind));
    
    ti = 'Average of Samples';
    FFdist_plots(ind+dg,mask_distance,maxlen,mon_ave,pol_ave,i_ave,y_ave,...
        std_ave,ti,0,seg_results);
     
end





%%


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
