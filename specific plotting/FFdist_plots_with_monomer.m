function FFdist_plots_with_monomer(ind,mask_distance,maxlen,mon,pol,intensity,y,stdpr,...
    ti,show_spindle,seg_results)

if ~show_spindle
    figure(ind);clf;
    %Plot monomer vs distance
    subplot(1,3,2);
    plot(mask_distance,mon,'go');
    xlabel('distance (microns)');
    ylabel('Monomer Concentration (au)');
    axis([min(mask_distance) maxlen 0 1]);
    %Plot polymer vs distance
    subplot(1,3,3);
    plot(mask_distance,pol,'ko');
    xlabel('distance (microns)');
    ylabel('Polymer Concentration (au)');
    axis([min(mask_distance) maxlen 0 1]);
    %Plot FRET fraction and intensity vs distance
    subplot(1,3,1);
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
    
    figure(ind);clf;
    %Plot monomer vs distance
    subplot(1,4,2);
    plot(mask_distance,mon,'go');
    xlabel('distance (microns)');
    ylabel('Monomer Concentration (au)');
    axis([min(mask_distance) maxlen 0 1]);
    %Plot polymer vs distance
    subplot(1,4,3);
    plot(mask_distance,pol,'ko');
    xlabel('distance (microns)');
    ylabel('Polymer Concentration (au)');
    axis([min(mask_distance) maxlen 0 1]);
    %Plot FRET fraction and intensity vs distance
    subplot(1,4,1);
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
    
    %seg_results = varargin{1};
    subplot(1,4,4); 
    pos_dis = find(seg_results.mask_distance>0);
    B= imoverlay(mat2gray(seg_results.image_stack),...
        seg_results.int_masks(:,:,pos_dis(1)));    
    imshow(B,'InitialMagnification','fit');
    title('Intensity image stack');
end


end