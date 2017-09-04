%%Used to plot polymer FF vs pixel group for spindle area FRET fraction
%%measurements.

clear;
zoom8x_April13 = [30635,30637];

% %bilinear and issue with 0-distance mask size
% zoom8x_Mar7 = [28519:28521,28523,28524];%bilinear
% zoom8x_Mar2 = 28510:28512;%this bilinear
% zoom8x_Jan26 = [28490:28493];%nearest [28498:28501];%bilinear: [28490,28491:28493];
% zoom12x = 28515:28517;%bilinear
% zoom8x = [zoom8x_Mar7,zoom8x_Mar2,zoom8x_Jan26];
% zoom2x = 28507:28508;%28506:28509;%bilinear

zoom12x_donor = 28518;
zoom8x_donor = [28522,28513:28514,28487,28527];%28527 is from 2-28, where there was no positive control
zoom2x_donor = [28525,28526];

zoom8x_Mar7 = [30720,30721,30723,30724];
zoom8x_Mar2 = 30726:30728;
zoom8x_Jan26 = 30731:30734;
zoom8x = [zoom8x_Mar7,zoom8x_Mar2,zoom8x_Jan26];
zoom12x = 30740:30742;
zoom2x = 30736:30739;

zoom8x_donor = [30735,30729,30730];%,30725];
zoom12x_donor = 30743;

%quad_div = 29449;
quad_div_1p5um =[30699:30702,30756:30763];
quad_div_1p5um = [quad_div_1p5um(1:9),quad_div_1p5um(11:12)];
quad_div = [30703:30706,30709:30711,30714:30718];

nospindle = [30648:30651,30653]; %too many pelicules
zoom8x_Apr13 = [30635,30637];
ivec = [28498:28501];%zoom8x_Jan26;%28330:28332; %
average_spindles = 1;
show_spindle = 0;
show_masks = 0;
showmon = 1;
pf = .1;
donor_offset = 1;



%hiint_28490
%hiint_28512
%[fret, fret_var, int, pol,mon,mdist] = FF_vs_dist_int([zoom8x_Mar2,zoom8x_Jan26],average_spindles,...
%    show_spindle,pf,showmon);
[fret, fret_var, int, pol,mon,mdist] = FF_vs_dist3(quad_div_1p5um,...
    average_spindles,show_spindle,pf,showmon,show_masks,donor_offset);

%%

%PLOT MON AND POL
% polave = mean(pol,1);
% polave_std_err = std(pol,1)/sqrt(11);
% monave = mean(mon,1);
% monave_std_err = std(mon,1)/sqrt(11);
% 
% msize_dot = 30;
% figure(1); clf; hold on;
% yyaxis left;
% errorbar(mdist,polave,polave_std_err,'.-','MarkerSize',msize_dot,'Color','k');
% ylabel('MT concentration (uM)');
% axis([-5 15 0 65]);
% yyaxis right;
% errorbar(mdist,monave,monave_std_err,'.--','MarkerSize',10,'Color',[.5,.5,.5]);
% ylabel('monomer concentration (uM)');
% xlabel('distance');
% title('Scaled by setting monomer to 18uM far from the spindle');
% axis([-5 15 0 65]);
% 
% 
% figure(2); clf; 
% 
% errorbar(mdist,polave,polave_std_err,'.-','MarkerSize',msize_dot,'Color','k');
% hold on; axis([-5 15 0 65]);
% errorbar(mdist,monave,monave_std_err,'.--','MarkerSize',10,'Color',[.5,.5,.5]);
% axis square

% 
%   csvwrite('distances',mdist);
%   csvwrite('polave',polave);
%   csvwrite('polave_std_err',polave_std_err);
%   csvwrite('monave',monave);
%   csvwrite('monave_std_err',monave_std_err);


%load(['/Users/bryankaye/Documents/MATLAB/data/matin/matin',num2str(i),'.mat']);
%imshow(mat2gray(seg_results.image_stack));

% grouping = [1,1,1,2,2,3,3,3,4,5,5,5];
% csvwrite('fret',fret);
% csvwrite('fret_var',fret_var);
% csvwrite('intensity',int);
% csvwrite('distances',mdist);
% csvwrite('group_key',grouping);

%[fret_donor,int_donor,pol,mon,mdist_donor] = FF_vs_dist(zoom8x_donor,average_spindles,show_spindle);



