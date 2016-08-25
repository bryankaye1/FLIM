clear;

jmax = 60;
base_filename = '2UF';
jind = 0;
combine_images =0;

for j = 2:jmax
    jind = jind +1;
    pth_data{jind} = 'C:\Users\Bryan\Documents\MATLAB\data\2016-07-31\';
    
    if jmax < 10
        dataname{jind} = strcat(base_filename,'_c',num2str(j));
    else
        if j < 10
            dataname{jind} = strcat(base_filename,'_c0',num2str(j));
        else
            dataname{jind} = strcat(base_filename,'_c',num2str(j));
        end
    end
end

[int_images] = sdt2image(pth_data, dataname);

if combine_images
    
    for j = 1:length(int_images)
        if j ==1
            mulexpim = int_images{j};
        else
            mulexpim = int_images{j} + mulexpim;
        end
    end
    %%
    blackval = min(min(mulexpim));
    whiteval = max(max(mulexpim));
    mulexpim = mat2gray(mulexpim,[blackval whiteval]);
    figure(21); clf; imshow(mulexpim); title('combined image')
end