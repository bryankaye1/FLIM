function [data_smoothed] = boxcar_averager_FLIM(data,reach)
data_smoothed = zeros(size(data));
[~,rows,cols] = size(data(1,:,:));
for i=1+reach:rows-reach
    for j = 1+reach:cols-reach
        for ip = -reach:reach
            for jp = -reach:reach                
                data_smoothed(:,i,j) = data(:,i+ip,j+jp)+data_smoothed(:,i,j);                
            end
        end
    end
end
end