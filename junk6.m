data = ones(5,6);
data_smoothed = zeros(size(data));
[rows,cols] = size(data);
reach = 1;
for i=1+reach:rows-reach
    for j = 1+reach:cols-reach
        for ip = -reach:reach
            for jp = -reach:reach                
                data_smoothed(i,j) = data(i+ip,j+jp)+data_smoothed(i,j);                
            end
        end
    end
end

data_smoothed