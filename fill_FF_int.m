function [x_int,x,y,stdpr,pixel_counts] = fill_FF_int(x_int,x,y,stdpr,output)

for fill_ind=output(1,1,1).jmax+1:output(1,1,1).num_int_bins
    x_int(fill_ind) = 0;
    x(fill_ind) = 0;
    y(fill_ind) = 0;
    stdpr(fill_ind) = 1e-10;
end

end