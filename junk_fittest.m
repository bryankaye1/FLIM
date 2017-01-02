

function [f] =  junk_fittest(p0,x,data)
      b = 0;%ones(length(x),1);
ft = fittype( 'bryansfn2(m,b,x)', 'problem', 'b','independent', 'x');

f = fit(x', data', ft, 'StartPoint', p0,'problem',b);


 
% ft = fittype( 'bryansfn2( m , x )' );
% f = fit( x', data', ft, 'StartPoint', 1);
end

