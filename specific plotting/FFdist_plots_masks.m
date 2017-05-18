function FFdist_plots_masks(ind,mask_distance,maxlen,pol,intensity,y,stdpr,...
    ti,show_spindle,seg_results)


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
    plot(mask_distance,pol,'ko');
    xlabel('distance (microns)');
    ylabel('Polymer Concentration (au)');
    axis([min(mask_distance) maxlen 0 1]);
    
    %Plot FRET fraction and intensity vs distance
    subplot(1,2,1);
    errorbar(mask_distance,y,stdpr,'bo');
    xlabel('distance (microns)');
    ylabel('FRET fraction');
    title(ti);
    axis([min(mask_distance) maxlen 0 .17]);
    % axis([-2 110 0 .13]);
    hold on;
    
    fret_peak =max(y);
    int_peak = max(intensity);
    intensity_norm = intensity * fret_peak/int_peak;
    plot(mask_distance,intensity_norm,'ro');
    legend('FRET-fraction','Intensity');
    
else
    %%
    color1 = [254,210,175]/255;
    color2 = [253,141,60]/255;
    color3 = [166,54,3]/255;
    
    green = [27,158,119]/255;%[217,95,2]/255;
    purple = [117,112,179]/255;
    pos_dis = find(seg_results.mask_distance>0);   
    msize_dot = 20;
    msize_aster = 10;
    
    figure(ind);clf;
%     %Plot monomer vs distance
%     subplot(1,4,2);
%     plot(mask_distance,mon,'go');
%     xlabel('distance (microns)');
%     ylabel('Monomer Concentration (au)');
%     axis([min(mask_distance) maxlen 0 1]);
    
    %Plot polymer vs distance
    subplot(3,1,3);
    yvec = 0:.01:2;
    plot(mask_distance,pol,'k.','MarkerSize',msize_dot);
    xlabel('distance (microns)');
    ylabel('Polymer Concentration (au)');
    axis([min(mask_distance)-0.5 maxlen+0.5 0 1.05]);
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
    axis([min(mask_distance)-0.5 maxlen+0.5 0 .16]);
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
    %seg_results = varargin{1};
   
    subplot(3,1,1); 
    pos_dis = find(seg_results.mask_distance>0);
    
    B= imoverlay(mat2gray(seg_results.image_stack),...
        seg_results.int_masks(:,:,pos_dis(1)-7),color3);
    
    B= imoverlay(B,seg_results.int_masks(:,:,pos_dis(1)-1),color2);
    
    B= imoverlay(B,seg_results.int_masks(:,:,pos_dis(1)+20),color1);
    
    imshow(B,'InitialMagnification','fit');
    title('Intensity image stack');
end


end