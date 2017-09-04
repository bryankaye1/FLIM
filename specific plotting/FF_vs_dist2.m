%%Used to plot polymer FF vs pixel group for spindle area FRET fraction
%%measurements.

function [ymat,varmat,imat,polmat,monmat,mask_distance] = FF_vs_dist2(ivec,...
    average_spindles,show_spindle,pf,showmon,show_masks,donor_offset)

green = [27,158,119]/255;%[217,95,2]/255;
purple = [117,112,179]/255;
quad_colors = {green,purple-.3,purple,green-[.1,.3,.3]};
quad_colors_std = {[194,213,204],[222,221,237],[204,202,218],[200,232,222]};
tub_con = 18;
ax_pol_hi = 80;

ind = 0;
load('/Users/bryankaye/Documents/MATLAB/FLIM/pf_est_table');
load('donor_offset.mat');
for i = ivec
    %%Load in FRET/intensity/distances
    
    clear x y distance stdpr
    ind = ind+1;
    [~,output] = load_mat_data(i,'local',1);
    load(['/Users/bryankaye/Documents/MATLAB/data/matin/matin',...
        num2str(i),'.mat']);
    if seg_results.input_params.scan_mag == 8
        maxlen = 17;%17.5;%17;         %15%20;
        minlen = -7;%-7.1384;%-5.5;        %-5.75;
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
    %ti = strrep([output(1,1,1).pth_data(end-10:end),output(1,1,1).dataname],'_',' ');
    ti = [];
    
    %%
    
    switch isfield(seg_results,'div_groups')
        case 0
            intensity = intensity(des_ind);
            y = y(des_ind);
            y0 = mean(FRET_offset);
            stdpr = stdpr(des_ind);

            des_row = find(ismember(pf_est.filename,...
                       [output(1,1,1).pth_data,output(1,1,1).dataname]));
            pf = pf_est.pf_fit(des_row)-y0;
            
            mon_raw = intensity.*(pf-y) ./ (pf*(1+(al-1).*y));
            eNmon = intensity.*(pf-(y-y0)) ./ (pf*(1+(al-1).*(y-y0)));
            %mon = mon/max(mon);
            
            if donor_offset
                [eNpol,polstd] = calc_pol(al,pf,y-mean(FRET_offset),stdpr,intensity);
            else
                [eNpol,polstd] = calc_pol(al,pf,y,stdpr,intensity);
            end
            ep = (eNmon + eNpol)/tub_con;
            pol_norm = eNpol/ ep;
            polstd_norm = polstd / ep;
            
            if average_spindles
                imat(ind,:) = intensity;
                ymat(ind,:) = y;
                polmat(ind,:) = pol_norm;
                monmat(ind,:) = eNmon;
                varmat(ind,:) = stdpr.*stdpr;
            else
                ymat = 0;
                varmat = 0;
                imat = 0;
                polmat = 0;
                monmat = 0;
            end
            
            if showmon
                %Plot monomer vs distance
                figure(ind);
                plot(mask_distance,mon,'bo');
                xlabel('distance (microns)');
                ylabel('Monomer Concentration (au)');
                axis([min(mask_distance) maxlen 0 1]);
                title(ti);
            end
            if ivec(end)==30731
                FFdist_plots(ind,mask_distance,maxlen,pol_norm,polstd_norm,...
                    intensity,y,stdpr,ti,show_spindle,seg_results,show_masks,0);
            end
            
            %     num_pho = intensity.*sum(sum(seg_results.int_masks));
            %     if sum(num_pho<500)
            %         fprintf('Photons count less than 500!'); pause;
            %     end
            
            if average_spindles
                imat(ind,:) = intensity;
                ymat(ind,:) = y;
                polmat(ind,:) = pol_norm;
                monmat(ind,:) = mon;
                varmat(ind,:) = stdpr.*stdpr;
                mask_distancemat(ind,:) = mask_distance;
            else
                ymat = 0;
                varmat = 0;
                imat = 0;
                polmat = 0;
                monmat = 0;
            end
            
        case 1 
            dgmax = length(seg_results.div_groups);
            for dg=1:dgmax
                dg_ind = des_ind+(dg-1)*length(seg_results.div_groups{1});
                intensity_dg = intensity(dg_ind);
                y_dg = y(dg_ind);
                stdpr_dg = stdpr(dg_ind);
                
                eNmon = intensity_dg.*(pf-y_dg) ./ (pf*(1+(al-1).*y_dg));
                %mon = mon/max(mon);
                
                [eNpol,polstd] = calc_pol(al,pf,y_dg-mean(FRET_offset),stdpr_dg,intensity_dg);
                ep = (eNmon + eNpol)/tub_con;
                pol_norm = pol/ ep;
                polstd_norm = polstd / ep;
                
                if show_spindle
                    % FFdist_plots(ind-1+dg,mask_distance,maxlen,pol_norm,...
                    %    polstd_norm,intensity_dg,y_dg,stdpr_dg,...
                    %    ti,show_spindle,seg_results,0,0)
                    
                    figure(ind-1+dg); clf; subplot(1,2,1);
                    bw_quads = sum(seg_results.div_masks{dg},3);
                    color_quads{dg}(:,:,1) = bw_quads*quad_colors{dg}(1);
                    color_quads{dg}(:,:,2) = bw_quads*quad_colors{dg}(2);
                    color_quads{dg}(:,:,3) = bw_quads*quad_colors{dg}(3);
                    imshow(color_quads{dg},'InitialMagnification','fit');
                    
                    msize_dot = 20;
                    ax_lim_pol = [min(mask_distance)-0.5 maxlen+0.5 0 ax_pol_hi];
                    subplot(1,2,2);
                    errorbar(mask_distance,pol_norm,polstd_norm,'Color',...
                        quad_colors{dg},'MarkerSize',msize_dot);
                    xlabel('distance (microns)');
                    ylabel('Polymer Concentration (au)');
                    axis(ax_lim_pol);
                    axis square;
                    
                    %  B= imoverlay(mat2gray(seg_results.image_stack),...
                    %     seg_results.div_masks{dg}(:,:,3),green);
                    %imshow(B,'InitialMagnification','fit');
                end
                
                imat{dg}(ind,:) = intensity_dg;
                ymat{dg}(ind,:) = y_dg;
                polmat{dg}(ind,:) = pol_norm;
                monmat{dg}(ind,:) = eNmon;
                varmat{dg}(ind,:) = stdpr_dg.*stdpr_dg;
                polmat{dg}(ind,:) = pol_norm;
                polmatstd{dg}(ind,:) = polstd_norm;
                mask_distancemat(ind,:) = mask_distance;
                
            end
    end
end

%%MAKE PLOTS OF 
if isfield(seg_results,'div_groups')
    if length(ivec)==1
        figure(ind+1+dgmax); clf;
        spindle_quads = 0*color_quads{1};
        for dg = 1:dgmax
            spindle_quads = color_quads{dg}+spindle_quads;
            %plot(mask_distance,polmat{dg},'Color',quad_colors{dg},'MarkerSize',msize_dot);
            errorbar(mask_distance,polmat{dg},polmatstd{dg},'Color',...
                quad_colors{dg},'MarkerSize',msize_dot); hold on;
            xlabel('distance (microns)');
            ylabel('Polymer Concentration (au)');
            axis(ax_lim_pol);
            axis square;
            
        end
        figure(ind+dgmax); imshow(spindle_quads);
        
    else
        msize_dot = 20;
        %ax_lim_pol = [min(mask_distance)-0.5 maxlen+0.1 0 1.05];
        ax_lim_pol = [-5 16.5 0 ax_pol_hi];
        %figure(length(ivec)+2+dgmax); 
        figure;clf;
        for dg = 1:dgmax
            for ind = 1:length(ivec)
           % errorbar(mask_distancemat(ind,:),polmat{dg}(ind,:),...
            %    polmatstd{dg}(ind,:),'Color',quad_colors{dg},...
            %    'MarkerSize',msize_dot); hold on;
            plot(mask_distancemat(ind,:),polmat{dg}(ind,:),'Color',quad_colors{dg},...
                'MarkerSize',msize_dot); hold on;
            xlabel('distance (microns)');
            ylabel('Polymer Concentration (au)');
            axis(ax_lim_pol);
            axis square;
            
            %calculate hwhm and half width at 10% max
            [MaxVal,maxi] = max(polmat{dg}(ind,:));
            hmi = find(polmat{dg}(ind,:)<=MaxVal/2);
            hwhm = mask_distancemat(ind,hmi(1))-mask_distancemat(ind,maxi);
            hwhmat{dg}(ind) = hwhm;
            
            tmi = find(polmat{dg}(ind,:)<=MaxVal/10);
            hwtm = mask_distancemat(ind,tmi(1))-mask_distancemat(ind,maxi);
            hwtmat{dg}(ind) = hwtm;           
            end
            
        end

            [hd(1),pd(1)] = ttest(hwhmat{1},hwhmat{2});
            [hd(2),pd(2)] = ttest(hwhmat{1},hwhmat{3});
            [hd(3),pd(3)] = ttest(hwhmat{4},hwhmat{2});
            [hd(4),pd(4)] = ttest(hwhmat{4},hwhmat{3});
            
            [hs(1),ps(1)] = ttest(hwhmat{2},hwhmat{3});
            [hs(2),ps(2)] = ttest(hwhmat{1},hwhmat{4});

            [hd(1),ptd(1)] = ttest(hwtmat{1},hwtmat{2});
            [hd(2),ptd(2)] = ttest(hwtmat{1},hwtmat{3});
            [hd(3),ptd(3)] = ttest(hwtmat{4},hwtmat{2});
            [hd(4),ptd(4)] = ttest(hwtmat{4},hwtmat{3});
            
            [hs(1),pts(1)] = ttest(hwtmat{2},hwtmat{3});
            [hs(2),pts(2)] = ttest(hwtmat{1},hwtmat{4});
        
    end
end

if average_spindles && ~isfield(seg_results,'div_groups') && length(ivec)>1
    if ~(ivec(1)==30735)
        y_ave = mean(ymat,1);
        i_ave = mean(imat,1);
        mon_ave = mean(monmat,1);
        std_ave = sqrt(mean(varmat,1)/(ind));
        mask_distance_ave = mean(mask_distancemat,1);
        donor = 0;
        
        polave = mean(polmat,1);
        polave_std = std(polmat,1);
        polave_norm = polave/max(polave);
        polave_std_norm = polave_std / max(polave);
        
        figure(ind+1);
        errorbar(mask_distance_ave,polave_norm,polave_std_norm);
        
    else
        
        %y_ave1 = mean(ymat(2:3,:),1);
        %y_ave2 = mean(ymat(1,:),1);
        %y_ave = (y_ave1 + 2*y_ave2 ) / 3;
        y_ave = mean(ymat,1);
        
        i_ave1 = mean(imat(2:3,:),1);
        i_ave2 = mean(imat(1,:),1);
        i_ave = (i_ave1 + i_ave2 ) / 2;
        
        polave = mean(polmat,1);
        polave_std = std(polmat,1);
        polave_norm = polave/max(polave);
        polave_std_norm = polave_std / max(polave);
        
        mon_ave = mean(monmat,1);
        std_ave = sqrt(mean(varmat,1)/(ind));
        donor=1;
        FRET_offset = y_ave;
        save('donor_offset.mat','FRET_offset');
        
    end
    %     if donor_offset
    %         [polave,polave_std] = calc_pol(al,pf,y_ave-mean(FRET_offset),std_ave,i_ave);
    %     else
    %         [polave,polave_std] = calc_pol(al,pf,y_ave,std_ave,i_ave);
    %     end
    FFdist_plots(ind+1,mask_distance,maxlen,polave_norm,polave_std_norm,...
        i_ave,y_ave,std_ave,ti,0,seg_results,show_masks,donor);
    
    if showmon
        figure(ind+1); clf;
        plot(mask_distance,mon_ave,'bo');
        title('monomer average');
        xlabel('distance (microns)');
        ylabel('Monomer Concentration (au)');
        axis([min(mask_distance) maxlen 0 1]);
    end
end

if  average_spindles && isfield(seg_results,'div_groups') && ind>1
%    figure;
    mask_distance_ave = mean(mask_distancemat,1);
    for dg = [4,3,1,2]
        polave = mean(polmat{dg},1);
        polave_norm{dg} = polave;
        
        %polave_norm{dg} = polave/max(polave);
        %polave_norm_std{dg} = std(polmat{dg},1)/max(polave);
        
%         % Shaded region plot styles
%         x = mask_distance_ave;
%         y = polave_norm{dg};
%         stdev = polave_norm_std{dg};
%         % I create some yu and yl here, for the example
%         yu = y+stdev;
%         yl = y-stdev;
%         fill([x fliplr(x)], [yu fliplr(yl)], quad_colors{dg}, 'linestyle', 'none')
%         axis(ax_lim_pol);
%         axis square;
%         hold on;
    end
    figure;
    for dg = 1:dgmax
        plot(mask_distance_ave,polave_norm{dg},'Color',quad_colors{dg},...
            'MarkerSize',msize_dot); hold on;
        xlabel('distance (microns)');
        ylabel('Polymer Concentration (au)');
        axis(ax_lim_pol);
        axis square;
        
    end   
end

end






%%

% Shaded region plot styles
%         x = mask_distance;
%         y = polave_norm;
%         % I create some yu and yl here, for the example
%         yu = y+polave_std_norm;
%         yl = y-polave_std_norm;
%         fill([x fliplr(x)], [yu fliplr(yl)], quad_colors{dg}, 'linestyle', 'none')
%         hold all
%         plot(x,y)

%fprintf('%s is %3.0f\n',ti,i);

%       des_row = find(ismember(pf_est.filename,...
%       [output(1,1,1).pth_data,output(1,1,1).dataname]));
% pf = max([pf_est.pf_fit(des_row)+pf_est.pf_fit_std(des_row),...
%      pf_est.max_FRET(des_row)+pf_est.max_FRET_std(des_row)]);

% pf = max([pf_est.pf_fit(des_row),pf_est.max_FRET(des_row)]);
%pf = pf_est.max_FRET(des_row);

% pf = pf_est.max_FRET(des_row)+pf_est.max_FRET_std(des_row);



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
