function  [pol,ave_FRET,ave_FRETerr] = ff2pol(x,y,stdpr,al,pf)

pol = sum( x.*y./( pf *( 1 + (al-1).*y) ));
ave_FRET = sum(y.*x./(1-y+al*y))/sum(x./(1-y+al*y));

yw = (x./(1-y+al*y))/sum(x./(1-y+al*y));
ywsig = (sqrt(x)./(1-y+al*y))/sum(x./(1-y+al*y));

FRETerr = y.*yw.*sqrt( (ywsig./yw).^2 + (stdpr./y).^2 );
ave_FRETerr = sqrt( sum(FRETerr.^2));

end