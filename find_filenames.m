function [tsm,dataname_cell] = find_filenames(cpath,base_name,tsm)

name_list  = ls([cpath,base_name,'*']);
if ispc
    dlmwrite('C:\Users\Bryan\Documents\MATLAB\FLIM\datanames_preview.txt',...
        name_list,'delimiter','');
    winopen('C:\Users\Bryan\Documents\MATLAB\FLIM\datanames_preview.txt');
else
    prev_pth ='/Users/bryankaye/Documents/MATLAB/FLIM/datanames_preview.txt';
    dlmwrite(prev_pth,name_list,'delimiter','');
    %fileattrib(prev_pth,'+r');
    system(['open -a TextWrangler ' ...
        '/Users/bryankaye/Documents/MATLAB/FLIM/datanames_preview.txt']);
end

size_name_list = size(name_list);
junk = input('Check file name list then press enter if ok');
fprintf('Continuing...\n');


if ispc
    for dc_ind = 1:size_name_list(1)
        name=strrep(name_list(dc_ind,:),' ','');
        dataname_cell{dc_ind} = name;
    end
else   
    namelist_str = dir([cpath,base_name,'*']);
    for dc_ind = 1:length(namelist_str)
       dataname_cell{dc_ind} = namelist_str(dc_ind).name;
    end
end

if tsm>0
    dataname_num = {};
    for i = 1:length(dataname_cell)
        dataname_num{end+1} = {dataname_cell{i}, i};
    end
    [start_nums,end_nums,dataname_cell] = find_series(dataname_num);
    tsm = end_nums-start_nums+1;
end
end
