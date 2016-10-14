function [new_filename] = remove_sdt(filename)

if length(filename)>3 %fixes filename if sdt is appended to file name
    if strcmp(filename(end-3:end),'.sdt')
        filename=filename(1:end-4);
    end
end
new_filename = filename;

end