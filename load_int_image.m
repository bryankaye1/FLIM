function [int_image,seg_results,varargout] = load_int_image(i)

i = num2str(i);
lastwarn('warning_reset_by_load_int_image.m');
warning('off','MATLAB:load:variableNotFound');

if ispc
    nstr_input = ['Y:\Users\bkaye\cluster\matin\matin',i,'.mat'];
else
    nstr_input = ['/Users/bryankaye/Documents/MATLAB/data/matin/matin',i,'.mat'];
end
load(nstr_input,'-mat','int_image','seg_results');
varargout{1} = 0;

[msgstr,~] = lastwarn;
if strcmp(msgstr,'Variable ''int_image'' not found.')
    fprintf('int_image%s DNE\n',i);
    varargout{1} = 1;
    int_image = [];
end
if strcmp(msgstr,'Variable ''seg_results'' not found.')
    fprintf('seg_results%s DNE\n',i);
    varargout{1} = 1;
    seg_results = 0;
end
warning('on','MATLAB:load:variableNotFound');

end










%     if length(seg_results)>1
%         varargout{1} = seg_results;
%       %  varargout{1} = seg_results;
%     else
%         varargout{1} = mat2gray(int_image);
%     end
% else
%     varargout{1} = mat2gray(int_image);
% end