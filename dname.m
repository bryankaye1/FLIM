function [dataname] = dname(basename, tseries, findw,expt,exptmax)
dataname = basename;
if findw > 1
    dataname = strcat(basename,num2str(expt));
end
if tseries > 0
    if exptmax < 10
        dataname = strcat(basename,'_c',num2str(expt));
    elseif exptmax < 100
        if expt<10
            dataname = strcat(basename,'_c0',num2str(expt));
        else
            dataname = strcat(basename,'_c',num2str(expt));
        end
    end
end

if tseries || findw > 0
    fprintf('dataname is %s\n', dataname);
    
end