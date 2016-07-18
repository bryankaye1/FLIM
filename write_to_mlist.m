function [MatName,SimName] = write_to_mlist(varargin)

numvarargs = length(varargin);
optargs = {0};
optargs(1:numvarargs) = varargin;
[specify_matin] = optargs{:};

clear mlist % I think this line does nothing

load('Y:\Users\bkaye\cluster\matlist.mat','-mat','mlist');

if round(specify_matin)==0
index = length(mlist.name)+1;
mlist.name(index) = index;
mlist.read(index) = 0;
else
index = round(varargin{1});
mlist.name(index) = index;
mlist.read(index) = 0;  
end

save('Y:\Users\bkaye\cluster\matlist.mat','mlist');
MatName = strcat('Y:\Users\bkaye\cluster\matin\matin',num2str(index),'.mat');
SimName = strcat('Y:\Users\bkaye\cluster\sim\sim',num2str(index),'.mat');
end