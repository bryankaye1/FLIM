% This m-file reads in rotation mount angles and power at the objective and
% outputs the interpolated power at the objective as a function of rotation
% mount angle


inrot = [3000,10000,15000,20000,30000,40000];
inpower = [100, 94,85,70,37,11];


inrot = (inrot/1000)*(pi*180); %Convert rotmount to angle then rads

cos2 = @ (A,theta_0,x) A*cos(x-theta_0).^2;
fresult = fit(inrot',inpower',cos2,'StartPoint',[max(inpower),0],...
    'Weights',1,'Lower',[0,-90],'Upper',[100,90]);

fresult