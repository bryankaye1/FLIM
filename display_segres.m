function [image2display] = display_segres(imagedata,seg_results)

if iscell(seg_results)
    
    [rows,cols] = size(seg_results{1});
    image2display = zeros(rows,cols);
    for i = 1:length(seg_results)
        image2display(1+rows*(i-1):rows*i,1:cols) = seg_results{i};
    end
    
else
    
    image2display = imagedata;
    
end