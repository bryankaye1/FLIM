
jmax = 10;
base_filename = 'm5';
for j = 1:jmax
    pth_data{j} = 'C:\Users\Bryan\Documents\MATLAB\data\2015-4-2\';
    
    if jmax < 10
        dataname{j} = strcat(base_filename,'_c',num2str(j));
    else
        if j < 10
            dataname{j} = strcat(base_filename,'_c0',num2str(j));
        else
            dataname{j} = strcat(base_filename,'_c',num2str(j));
        end
    end
end

[~] = sdt2image(pth_data, dataname);