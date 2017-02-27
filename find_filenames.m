function [tsm,dataname_cell] = find_filenames(cpath,base_name,tsm,varargin)

%If varargin{1} exists, then the script will not ask you to confirm
%filenames are correct

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
if length(varargin)==0
junk = input('Check file name list then press enter if ok');
fprintf('Continuing...\n');
end


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
    base_filenames_newline = '';
    base_filenames_comma = '';
    for i = 1:length(dataname_cell)
        base_filenames_newline = sprintf('%s %s \n',base_filenames_newline,dataname_cell{i});
        base_filenames_comma = sprintf('%s''%s'',',base_filenames_comma,dataname_cell{i});
    end
    tsm = end_nums-start_nums+1;
    %For copying the filenames for pasting directly into code
    base_names_pth ='/Users/bryankaye/Documents/MATLAB/FLIM/combined_datanames.txt';
    dlmwrite(base_names_pth,[base_filenames_newline,base_filenames_comma],'delimiter','');
    system(['open -a TextWrangler ' base_names_pth]);
    bryan = 1;
    
    
    
end
end
