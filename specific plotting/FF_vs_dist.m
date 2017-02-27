%%Used to plot polymer FF vs pixel group for spindle area FRET fraction
%%measurements.
clear;
ivec = 28303:28306;
ind = 0;
average_spindles = 1;
maxlen = 0;
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

    %%Get FRET fraction (fluoro population) and intensity
  
    for j = 1:length(output)
        x(j) = ni(j)*(sum(output(1,1,j).prest'.*output(1,1,j).prestx)...
            +sum(output(1,1,j).w02est'.*output(1,1,j).w02estx)); %#ok<SAGROW>
        [y(j),stdpr(j)] = transform_wf_to_f(output(1,1,j).prest,...
            output(1,1,j).prestx,al,output(1,1,j).w02est,...
            output(1,1,j).w02estx,x(j),'dont_combine'); %#ok<*SAGROW>
    end  
    
    num_pixels = squeeze(sum(sum(int_masks,1),2))';
    intensity = x./num_pixels;   
    figure(ind); subplot(1,2,1);
    plot(mask_distance,y,'bo')
    %errorbar(distance,y,stdpr);
    xlabel('distance');
    ylabel('FRET fraction');
    title(strrep(output(1,1,1).dataname,'_',' '));
   % axis([-5 20 0 .18]);
   %axis([-10 110 0 .13]);
    
    
    subplot(1,2,2); %errorbar(distance,intensity,sqrt(intensity));
    plot(mask_distance,intensity,'ro');
    xlabel('distance');
    ylabel('intensity');
    title(strrep(output(1,1,1).dataname,'_',' '));
    
    if average_spindles
        xmat{ind} = mask_distance;
        ymat{ind} = y; 
        len = length(x);
        if len > maxlen
        maxlen = len; 
        max_len_ind = ind;
        end
    end 
end

x_sq = nan(ind,max(maxlen));
y_sq = nan(ind,max(maxlen));

for k = 1:ind
x_sq(k,1:length(xmat{k})) = xmat{k};
y_sq(k,1:length(ymat{k})) = ymat{k};
end

x_ave = nanmean(x_sq,1);
y_ave = nanmean(y_sq,1);

figure(ind+1); plot(x_ave,y_ave,'bo');
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
