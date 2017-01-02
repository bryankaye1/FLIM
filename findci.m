%This function starts at the peak of the PDF and moves outward in equal
%steps left and right, until the area under the curve is equal to "desci".
%If the left or right bound hits the edge of the PDF, the other bound keeps
%continuing outward until area under the curve is desci. This function will
%print an error to the screen if both bounds hit the end points of the
%function before the area under the curve reaches desci.

%Inputs: 
% pdf: y component of probability function (unnormalized ok). 
% pdfx: x compenent of probability function (has scale) 
% desci: (desired confidence interval) Area under curve corresponging to error.
%   If set to .68, then this equals one standard deviation (unnormalized NOT ok)
%
%varargin {1}: If set to 'confidence_bounds', the function will output the
%lower and upper limit to the confidence interval. 

%If set to 'error_size', it outputs the lower error size (distance form mode
%to leftmost confidence bound) and upper error size. Useful for use in 
%Matlab's 'error' plotting function.

%If set to 'relative_endpoints', it outputs the relative probability of the
%endpoints relative to the max of the pdf

%outputs: confidence interval (standard deviation), cimax, cimin


function [cint,varargout] = findci(pdf, pdfx, desci,varargin)


numvarargs = length(varargin);
if numvarargs > 1
    error('myfuns:findci:TooManyInputs','requires at most 1 optional inputs');
end

optargs = {'vargarin_default'};
optargs(1:numvarargs) = varargin;

pdf = pdf/sum(pdf); %Normalize PDF

[~,mode] = max(pdf); %find mode
aci = 0;
lc = 0;
rc =0;
stopr=0;
stopl = 0;

while aci < desci
    
    if mode + rc ==length(pdf)
        stopr = 1;
    else
        rc =rc+1;
    end
    
    if mode - lc == 1;
        stopl = 1;
    else
        lc = lc+1;
    end
    
    aci = sum(pdf(mode-lc:mode+rc));
    if stopl ==1 && stopr==1
        fprintf('error in cint function');
        break;
    end
    
end
cint = (pdfx(mode+rc)-pdfx(mode))/2;
if strcmp(optargs,'confidence_bounds')
    varargout{1} = pdfx(mode-lc);
    varargout{2} = pdfx(mode+rc);
elseif strcmp(optargs,'error_size')
    varargout{1} = pdfx(mode)-pdfx(mode-lc);
    varargout{2} = pdfx(mode+rc)-pdfx(mode);
elseif strcmp(optargs,'dist_endpoints')
    varargout{1} = pdf(mode)/pdf(1);
    varargout{2} = pdf(mode)/pdf(end);
    
end
end