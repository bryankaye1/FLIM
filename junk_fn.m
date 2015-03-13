function [out1,out2,varargout] = junk_fn(a,b,varargin)
if nargin > 2
out2 = varargin{2};
out1 = a + b + varargin{1};
else
   out2 = 'no 2nd input'; 
   out1 = a + b ;
end
if a+b ==pi
    varargout{1} = [1:10];
end

end