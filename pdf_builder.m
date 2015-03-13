function [irf_pdf] = pdf_builder(fun)
%T is the domain (range in X) of the function. For FLIM, T is the excitation
%time
%fun is the function
T = 12.58;
f = (10^7)*fun/sum(fun);

dist = [];
for i=1:length(f)
    
    t = i*(T/length(f))*ones(1,round(f(i)));
    dist=[dist t];
    
end
irf_pdf = dist;
end

% To draw from this distribution
% data = pdf_builder(f)
% randraws = data(randi(length(data),1,1000));