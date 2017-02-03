function  [pol,ave_FRET,ave_FRETerr,varargout] = ff2pol(x,y,stdpr,al,pf,varargin)

%vararg to compute fraction of donor in polymer by sum(pol) / sum(pol+mon)

pol = sum( x.*y./( pf *( 1 + (al-1).*y) ));
if length(varargin)
pol_fraction = pol / (pol+length(y)*varargin{1});
varargout{1} = pol_fraction;
end

ave_FRET = sum(y.*x./(1-y+al*y))/sum(x./(1-y+al*y));
%where did this come from? It looks like the above equation but without the
%divided by Pf.

yw = (x./(1-y+al*y))/sum(x./(1-y+al*y));
ywsig = (sqrt(x)./(1-y+al*y))/sum(x./(1-y+al*y));

FRETerr = y.*yw.*sqrt( (ywsig./yw).^2 + (stdpr./y).^2 );
ave_FRETerr = sqrt( sum(FRETerr.^2));



end