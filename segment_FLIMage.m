function [outFLIM,outint,out_image] = segment_FLIMage(inint,inFLIM,blur_dist,threshold)   

int_image = mat2gray(make2D(inint)); %make square on int_image? Can you imgaussfilt a line?
figure(1); subplot(2,2,1); imshow(int_image,'InitialMagnification', 'fit');

blurred = imgaussfilt(int_image,blur_dist);
subplot(2,2,2);
imshow(blurred,'InitialMagnification', 'fit');

thr1 = int_image;
thr1(thr1 >= threshold) = 1;
thr1(thr1 < threshold) = 0;
subplot(2,2,3);
imshow(thr1,'InitialMagnification', 'fit');

thr3 = blurred;
thr3(thr3 >= threshold) = 1;
thr3(thr3 < threshold) = 0;
subplot(2,2,4); imshow(thr3,'InitialMagnification', 'fit');
thr3l = logical(thr3);

outint = inint(thr3l);
outFLIM = inFLIM(thr3l,:);
out_image={int_image,blurred,thr1,thr3};
drawnow;
end

function [square_image] = make2D(oimage)
square_image = reshape(oimage,128,round(length(oimage)/128));
end

% function [line_image] = make1D(oimage)
% line_image = reshape(oimage,1,
% end


% for l = 1:128*length(inint_line)
%     j = mod(l-1,128)+1;
%     k = floor((l-1)/128) + 1;   
%     inint(j,k) = inint_line;
%     flimap(j,k)= fBest;
% end