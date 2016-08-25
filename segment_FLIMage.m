function [outint,outFLIM] = segment_FLIMage(inint,inFLIM,blur_dist,threshold)

int_image = mat2gray(inint);
figure(1); subplot(2,2,1); imshow(makesquare(int_image));

blurred = imgaussfilt(int_image,blur_dist);
subplot(2,2,2);
imshow(makesquare(oimage));

thr1 = int_image;
thr1(thr1 >= threshold) = 1;
thr1(thr1 < threshold) = 0;
subplot(2,2,3);
imshow(makesqaure(thr1));

thr3 = blurred;
thr3(thr3 >= threshold) = 1;
thr3(thr3 < threshold) = 0;
subplot(2,2,4); imshow(makesquare(thr3));

outint = inint(thr3);
outFLIM = inFLIM(thr3);

end

function [square_image] = makesquare(oimage)
square_image = reshape(oimage,sqrt(length(oimage)),sqrt(length(oimage)));
end