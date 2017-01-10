function [x_int,x] = npho2cps(x,total_pixels,acq_time,output)


if isfield(output,'pixel_counts')
    x_int = x.*(total_pixels/acq_time); % x = x.*(128*128/20); % x = ph/pixel  .* pixels/sec
    x = output(1,1,1).pixel_counts.*x; %pixel * ph/pixel
else
    delta_t = acq_time / length(output(1,1,1).ni);
    x_int = x./delta_t;
end

end