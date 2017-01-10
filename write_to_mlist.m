function [MatName,index] = write_to_mlist(varargin)

numvarargs = length(varargin);
optargs = {0};
optargs(1:numvarargs) = varargin;
[specify_matin] = optargs{:};

clear mlist % I think this line does nothing
if ispc
    load('Y:\Users\bkaye\cluster\matlist.mat','-mat','mlist');
else
    load('/Users/bryankaye/Dropbox/matlist.mat','-mat','mlist');
end

if round(specify_matin)==0
    index = length(mlist.name)+1;
    mlist.name(index) = index;
    mlist.read(index) = 0;
else
    index = round(varargin{1});
    mlist.name(index) = index;
    mlist.read(index) = 0;
end

if ispc
    save('Y:\Users\bkaye\cluster\matlist.mat','mlist');
    MatName = strcat('Y:\Users\bkaye\cluster\matin\matin',num2str(index),'.mat');
else
    save('/Users/bryankaye/Dropbox/matlist.mat','mlist');
    MatName = strcat('/Users/bryankaye/Documents/MATLAB/data/matin/matin',num2str(index),'.mat');

end
end