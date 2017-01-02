function[start_num,end_num,base_name] = find_series(name_num)
%remove .sdt from filenames
for i = 1:length(name_num)
   name_num{i}{1} = strrep(name_num{i}{1},'.sdt','');
end

%Find every the time-series and record their names
base_name = {};
base_name_ind = [];
start_num = [];
for i = 1:length(name_num)
    ind = regexp(name_num{i}{1},'_c01$', 'once');
    if isempty(ind)
      ind = regexp(name_num{i}{1},'_c1$', 'once');
    end
    if ~isempty(ind)
        base_name{end+1} = name_num{i}{1}(1:ind-1);
        base_name_ind(end+1) = ind;
        start_num(end+1) = name_num{i}{2};
    end
end

%Find the last time point in each time-series
end_num = [];
for j = 1:length(base_name) %loop of series names
    largest = 0;
    for i = 1:length(name_num) %check each name to jth series name
        ind = regexp(name_num{i}{1},[base_name{j},'_c'], 'once');
        if ~isempty(ind)
            time_point = str2num(name_num{i}{1}(base_name_ind(j)+2:end));
            if time_point>largest
                largest = time_point;
                end_num(j) = name_num{i}{2};
            end
        end
    end
end

% test name_num below:
% name_num{1} = {'base_c01',1};
% name_num{2} = {'base_c02',2};
% name_num{3} = {'base_c03',3};
% name_num{4} = {'base2_c01',6};
% name_num{5} = {'base2_c02',7};
% name_num{6} = {'base2_c03',8};
% name_num{7} = {'base3',9};
% name_num{8} = {'base30_c01',10};
% name_num{9} = {'base30_c02',11};