function [tsm,dataname_cell] = find_filenames(cpath,base_name,tsm)

name_list  = ls([cpath,base_name,'*']);
dlmwrite('C:\Users\Bryan\Documents\MATLAB\FLIM\datanames_preview.txt',name_list,'delimiter','');
winopen('C:\Users\Bryan\Documents\MATLAB\FLIM\datanames_preview.txt')
size_name_list = size(name_list);
junk = input('Check file name list then press enter if ok');
fprintf('Continuing...\n');
for dc_ind = 1:size_name_list(1)
    name=strrep(name_list(dc_ind,:),' ','');
    dataname_cell{dc_ind} = name;
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