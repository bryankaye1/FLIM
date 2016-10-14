function [file_name2] = append_timeseries_names(file_name,ts,tsm)

file_name = remove_sdt(file_name); %removes ".sdt" from end of filename if present
if tsm ==1
    file_name2 = strcat(file_name,'.sdt');
else
    if tsm <10
        file_name2 = strcat(file_name,'_c0',num2str(ts),'.sdt');
    elseif tsm < 100
        if ts<10
            file_name2 = strcat(file_name,'_c0',num2str(ts),'.sdt');
        else
            file_name2 = strcat(file_name,'_c',num2str(ts),'.sdt');
        end
    end
end
end