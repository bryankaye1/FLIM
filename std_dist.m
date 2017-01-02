function [stdev] = std_dist(px,x)
stdev = sqrt ( sum(px.*x.^2) - sum(px.*x).^2 );
end