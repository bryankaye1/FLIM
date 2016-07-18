function [image_cell] = sdt2image(paths,file_names,varargin)
%This function converts SDT files to intensity images. The white value is 
%set to the number of photons in the brightest pixel of all the images 
%passed to this function; the black value is set to the number of photons 
%in the dimmest pixel.
%
% Inputs: 
% (1) cell of file paths (string).
% (2) cell of file names (string).
% (3) Vararg 1: if set to 'nofig', this will plot the images (but it will
% still return images).

% Outputs:
% (1) cell of intensity matrices (float, cell of matrices)
% Creates figures for each intensity

%Notes on setting defaults of varargins: http://blogs.mathworks.com/loren/2009/05/05/nice-way-to-set-function-defaults/




% only want 1 optional input at most
numvarargs = length(varargin);
if numvarargs > 1
    error('myfuns:sdt2image:TooManyInputs','requires at most 1 optional inputs');
end

optargs = {'n/a'};
optargs(1:numvarargs) = varargin;
make_figures = optargs;    

for j= 1:length(file_names)
    
    [pmat,jmax,ni,imagedata] = spc_2_his(1,4096,file_names{j},paths{j},1,1,'imout');
    raw_images{j} = squeeze(imagedata);
end

whiteval = 0;
blackval = inf;
for j = 1:length(file_names)
    whiteval_temp = max(max(raw_images{j}));
    blackval_temp = min(min(raw_images{j}));
    
    if  whiteval_temp > whiteval
        whiteval = whiteval_temp;
    end
    
    if  blackval_temp < blackval
        blackval = blackval_temp;
    end
end
%%
set(0,'DefaulttextInterpreter','none');
iptsetpref('ImshowInitialMagnification', 'fit');
for j = 1:length(file_names)
    image_cell{j} = mat2gray(raw_images{j},[blackval whiteval]);
    if strcmp(make_figures, 'nofig')
    else
    figure(j); clf; imshow(image_cell{j}); title(file_names{j});
    end
end

end