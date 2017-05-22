function FFdist_plots(ind,mask_distance,maxlen,pol,polstd,intensity,y,stdpr,...
    ti,show_spindle,seg_results,show_masks,donor)

color1 = [254,210,175]/255;
color2 = [253,141,60]/255;
color3 = [166,54,3]/255;

green = [27,158,119]/255;%[217,95,2]/255;
purple = [117,112,179]/255;
pos_dis = find(seg_results.mask_distance>0);
msize_dot = 20;
msize_aster = 10;
yvec = 0:.01:2;
ax_lim_FRET = [min(mask_distance)-0.5 maxlen+0.5 0 .16];
ax_lim_pol = [min(mask_distance)-0.5 maxlen+0.5 0 1.05];

if ~show_spindle
    figure(ind);clf;
    %Plot monomer vs distance
%     subplot(1,3,2);
%     plot(mask_distance,mon,'go');
%     xlabel('distance (microns)');
%     ylabel('Monomer Concentration (au)');
%     axis([min(mask_distance) maxlen 0 1]);

    %Plot polymer vs distance
    subplot(1,2,2);
    errorbar(mask_distance,pol,polstd,'k.','MarkerSize',msize_dot);
    xlabel('distance (microns)');
    ylabel('Polymer Concentration (au)');
    axis(ax_lim_pol);
    axis square;
    
    
    %Plot FRET fraction and intensity vs distance
    subplot(1,2,1);     
    errorbar(mask_distance,y,stdpr,'.','MarkerSize',msize_dot,'Color',purple);
    hold on;
    if donor
    fret_peak =0.1267;%*(549.5909/370.9563); %max FRET value from matout 30731
    int_peak = 549.5909;  
    else
    fret_peak =max(y);
    int_peak = max(intensity);
    end

    intensity_norm = intensity * fret_peak/int_peak;
    plot(mask_distance,intensity_norm,'*','MarkerSize',msize_aster,'Color',green);
    
   % legend('FRET-fraction','Intensity');
    title(ti);
    xlabel('distance (microns)');
    ylabel('FRET fraction');
    axis(ax_lim_FRET);
    axis square;


    
elseif show_masks  
    figure(ind);clf;
    
    %Plot polymer vs distance
    subplot(3,1,3);
    yvec = 0:.01:2;
    errorbar(mask_distance,pol,polstd,'k.','MarkerSize',msize_dot);
    xlabel('distance (microns)');
    ylabel('Polymer Concentration (au)');
    axis(ax_lim_pol);
    hold on;
    plot(mask_distance(pos_dis(1)-7)*ones(1,length(yvec)),yvec,'--',...
        'Color', color3,'LineWidth',2);
    plot(mask_distance(pos_dis(1)-1)*ones(1,length(yvec)),yvec,'--',...
        'Color', color2,'LineWidth',2);
    plot(mask_distance(pos_dis(1)+20)*ones(1,length(yvec)),yvec,'--',...
        'Color', color1,'LineWidth',2);
    axis square;
    
    %Plot FRET fraction and intensity vs distance
    subplot(3,1,2);
    errorbar(mask_distance,y,stdpr,'.','MarkerSize',msize_dot,'Color',purple);
    xlabel('distance (microns)');
    ylabel('FRET fraction');
    title(ti);
    axis(ax_lim_FRET);
    % axis([-2 110 0 .13]);
    hold on;
    
    fret_peak =max(y);
    int_peak = max(intensity);
    intensity_norm = intensity * fret_peak/int_peak;
    plot(mask_distance,intensity_norm,'*','MarkerSize',msize_aster,'Color',green);

    
    plot(mask_distance(pos_dis(1)-7)*ones(1,length(yvec)),yvec,'--',...
        'Color', color3,'LineWidth',2);
    plot(mask_distance(pos_dis(1)-1)*ones(1,length(yvec)),yvec,'--',...
        'Color', color2,'LineWidth',2);
    plot(mask_distance(pos_dis(1)+20)*ones(1,length(yvec)),yvec,'--',...
        'Color', color1,'LineWidth',2);
    legend('FRET-fraction','Intensity');
    axis square;
   
    subplot(3,1,1); 
    B= imoverlay(mat2gray(seg_results.image_stack),...
        seg_results.int_masks(:,:,pos_dis(1)-7),color3);
    
    B= imoverlay(B,seg_results.int_masks(:,:,pos_dis(1)-1),color2);
    
    B= imoverlay(B,seg_results.int_masks(:,:,pos_dis(1)+20),color1);
    
    imshow(B,'InitialMagnification','fit');
    title('Intensity image stack');
else
    
    
    figure(ind);clf;
    %     %Plot monomer vs distance
    %     subplot(1,4,2);
    %     plot(mask_distance,mon,'go');
    %     xlabel('distance (microns)');
    %     ylabel('Monomer Concentration (au)');
    %     axis([min(mask_distance) maxlen 0 1]);
    
   subplot(1,3,3);
    errorbar(mask_distance,pol,polstd,'k.','MarkerSize',msize_dot);
    xlabel('distance (microns)');
    ylabel('Polymer Concentration (au)');
    axis(ax_lim_pol);
    axis square;
    %Plot FRET fraction and intensity vs distance
    subplot(1,3,2);
    errorbar(mask_distance,y,stdpr,'.','MarkerSize',msize_dot,'Color',purple);
    xlabel('distance (microns)');
    ylabel('FRET fraction');
    title(ti);
    axis(ax_lim_FRET);
    % axis([-2 110 0 .13]);
    hold on;
    
    fret_peak =max(y);
    int_peak = max(intensity);
    intensity_norm = intensity * fret_peak/int_peak;
    plot(mask_distance,intensity_norm,'*','MarkerSize',msize_aster,'Color',green);
    legend('FRET-fraction','Intensity');
        axis square;
    
    %seg_results = varargin{1};
    subplot(1,3,1); 
    pos_dis = find(seg_results.mask_distance>0);
    B= imoverlay(mat2gray(seg_results.image_stack),...
        seg_results.int_masks(:,:,pos_dis(1)));    
    imshow(B,'InitialMagnification','fit');
    title('Intensity image stack');
end


end