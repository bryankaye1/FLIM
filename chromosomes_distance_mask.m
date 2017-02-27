
close all;
clear

data_pth = 'C:\Users\BIGSAS\Desktop\FLIM\';
data_name ='20170211_test.mat';
filename=[data_pth,data_name];
load(filename,'ellipse','s');

%%% CREATE DISTANCE MASK FROM CHROMOSOME BASED ON INNER-ELLIPSE
for k=1:size(ellipse,2)
    
    mask_dist{k}=im2double(ellipse{k}(:,:,1));
    mask_dist{k}(mask_dist{k}==1)=NaN;
    for i = (65-round(s{k}.MinorAxisLength/2)):1:(65+round(s{k}.MinorAxisLength/2))
       %for j=1:size(int_final{k},2)
       for j=1:size(ellipse{k},2)
            if ~isnan(mask_dist{k}(i,j))
                mask_dist{k}(i,j)=abs(j-65);
            end    
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     for i = 1:(65-round(s{k}.MinorAxisLength/2));
       for j=1:size(ellipse{k},2)
            if ~isnan(mask_dist{k}(i,j))
                mask_dist{k}(i,j)=sqrt((i-(65-s{k}.MinorAxisLength/2))^2+(j-65)^2);
           end    
        end
     end
    
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
   
   for i = (65+round(s{k}.MinorAxisLength/2)):1:size(ellipse{k},2);
       for j=1:size(ellipse{k},2)
            if ~isnan(mask_dist{k}(i,j))
                mask_dist{k}(i,j)=sqrt((i-(65+s{k}.MinorAxisLength/2))^2+(j-65)^2);
           end    
        end
   end
 
figure
imagesc(mask_dist{k});    %%%% this is the mask of equal distances
   
end



%%%% SORT DISTANCES FROM MAP TO LISTS
for k=1:size(ellipse,2)
for m=1:floor(max(max(mask_dist{k})))
    b=[];
    for i=1:size(ellipse{k},1)
        for j=1:size(ellipse{k},2)
            
            if m == floor(mask_dist{k}(i,j));
                
                coordinates=[j,i];
                b=[coordinates; b];
                 end    
        end  
    end
    
  list{k,m}=b;  
end
end
%%% list contains coordinates of pixels at same distances - due to
%%% non-integer distances - 1 pixel refers to the distance interval 1-1.99
%%% and so on; 1st coordiate refers to from left to right; 2nd from top to
%%% botoom;