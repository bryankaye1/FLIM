function [input] = ini_input(cyclesmax,exptmax,jmax)

input(cyclesmax,exptmax,jmax).jmax = pi;
input(cyclesmax,exptmax,jmax).exptmax = pi;
input(cyclesmax,exptmax,jmax).cyclesmax = pi;

input(cyclesmax,exptmax,jmax).datahis = pi*ones(1,2928);
input(cyclesmax,exptmax,jmax).ga = pi*ones(1,3681); %ga is name of vector of shifted irf

input(cyclesmax,exptmax,jmax).w1step = pi; input(cyclesmax,exptmax,jmax).w1min = pi; input(cyclesmax,exptmax,jmax).w1max = pi;
input(cyclesmax,exptmax,jmax).w2step = pi; input(cyclesmax,exptmax,jmax).w2min = pi; input(cyclesmax,exptmax,jmax).w2max = pi;
input(cyclesmax,exptmax,jmax).prstep = pi; input(cyclesmax,exptmax,jmax).prmin = pi; input(cyclesmax,exptmax,jmax).prmax = pi;
input(cyclesmax,exptmax,jmax).w02step = pi; input(cyclesmax,exptmax,jmax).w02min = pi; input(cyclesmax,exptmax,jmax).w02max = pi;
input(cyclesmax,exptmax,jmax).extstep = pi; input(cyclesmax,exptmax,jmax).extmin = pi; input(cyclesmax,exptmax,jmax).extmax = pi;
input(cyclesmax,exptmax,jmax).fracstep = pi;

input(cyclesmax,exptmax,jmax).dataname = 'initialized_name';
input(cyclesmax,exptmax,jmax).pth_data = 'initialized_path';
input(cyclesmax,exptmax,jmax).irf_name= 'initialized_name';
input(cyclesmax,exptmax,jmax).pth_irf = 'initialized_path';

input(cyclesmax,exptmax,jmax).data_shift_name = 'initialized_name';
input(cyclesmax,exptmax,jmax).pth_data_for_shift = 'initialized_name';


input(cyclesmax,exptmax,jmax).thrt = pi; %0 means do not threshold data, 1 means threshold , 2 means treshold for desired photon amount
input(cyclesmax,exptmax,jmax).imthr = pi; %if input.thresh/thrt is 1, thresh is amount to threshol by. Thresh is max(int)/imthr
input(cyclesmax,exptmax,jmax).DesPho = pi;  %if input.thresh/thrt is 2, DesPho sets the number of photons

input(cyclesmax,exptmax,jmax).thr = pi;
input(cyclesmax,exptmax,jmax).brem = pi;
input(cyclesmax,exptmax,jmax) .bins= pi;
input(cyclesmax,exptmax,jmax).tmini = pi;
input(cyclesmax,exptmax,jmax).tmaxi = pi;
input(cyclesmax,exptmax,jmax).ext= pi;
input(cyclesmax,exptmax,jmax).wig = pi;

input(cyclesmax,exptmax,jmax).sfracv = pi; input(cyclesmax,exptmax,jmax).shiftmin = pi; input(cyclesmax,exptmax,jmax).shiftmax= pi;
input(cyclesmax,exptmax,jmax).w2step_shift = pi; input(cyclesmax,exptmax,jmax).w2min_shift = pi; input(cyclesmax,exptmax,jmax).w2max_shift = pi;
input(cyclesmax,exptmax,jmax).backstep = pi; input(cyclesmax,exptmax,jmax).backmin = pi; input(cyclesmax,exptmax,jmax).backmax = pi;
input(cyclesmax,exptmax,jmax).wigstep = pi; input(cyclesmax,exptmax,jmax).wigmin = pi; input(cyclesmax,exptmax,jmax).wigmax = pi;

input(cyclesmax,exptmax,jmax).pth_wigs = 'initialized_path';
input(cyclesmax,exptmax,jmax).wigsname = 'initialized_name';
input(cyclesmax,exptmax,jmax).pth_ext = 'initialized_path';
input(cyclesmax,exptmax,jmax).extname = 'initialized_name';
end