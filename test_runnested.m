

function [f] =  junk_fittest(p0,x,data)
      b = 0;
      c = 1;
ft = fittype( 'bryansfn(  m , b, c, x)' , 'problem', {'b','c'}, 'independent','x');

f = fit(x', data', ft, 'StartPoint', p0,'problem',{'b','c'});


% x = 1:10;
% y = rand(1,10)+x;
% 
% ft = fittype( 'bryansfn( m , x )' );
% f = fit( x', y', ft, 'StartPoint', 1);
end

