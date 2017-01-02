function [input]  = contruct_input_struct(split_matin,jmax,exptmax,...
    cyclesmax,p,gab,w1step,w1min,w1max,w2step,w2min,w2max,prstep,prmin,prmax,...
    w02step,w02min,w02max,fracstep,shift,shiftb,dataname,pth_data,irfname,pth_irf,...
    data_shift_name,pth_data_for_shift,ngr,ni,thr,bneed,pulsewb,tmini,tmaxi,...
    ext,wigsb,pth_wigs,wigsname,pth_ext,extname,comment,tbac,tfw,reach,...
    combine_exposures,tsm,cindex,expt,jind,input)

if split_matin==0
    split_matin=1;
end

input(cindex,expt,jind).datahis = p;
input(cindex,expt,jind).ga = gab; %ga is name of vector of shifted irf

input(cindex,expt,jind).w1step = w1step;
input(cindex,expt,jind).w1min = w1min; input(cindex,expt,jind).w1max = w1max;
input(cindex,expt,jind).w2step = w2step;
input(cindex,expt,jind).w2min = w2min; input(cindex,expt,jind).w2max = w2max;
input(cindex,expt,jind).prstep = prstep;
input(cindex,expt,jind).prmin = prmin; input(cindex,expt,jind).prmax = prmax;
input(cindex,expt,jind).w02step = w02step;
input(cindex,expt,jind).w02min = w02min; input(cindex,expt,jind).w02max = w02max;
input(cindex,expt,jind).extstep = fracstep;
input(cindex,expt,jind).extmin = 0; input(cindex,expt,jind).extmax = 0;
input(cindex,expt,jind).fracstep = fracstep;

input(cindex,expt,jind).dataname = dataname;
input(cindex,expt,jind).pth_data = pth_data;
input(cindex,expt,jind).irf_name= irfname;
input(cindex,expt,jind).pth_irf = pth_irf;

if split_matin>1 || (cindex==1 && (expt==1 && jind==1))
    
    input(cindex,expt,jind).jmax = jmax/split_matin;
    input(cindex,expt,jind).exptmax = exptmax;
    input(cindex,expt,jind).cyclesmax = cyclesmax;
    
    input(cindex,expt,jind).backstep = shift.backstep;
    input(cindex,expt,jind).backmin = shift.backmin;
    input(cindex,expt,jind).backmax = shift.backmax;
    input(cindex,expt,jind).w2step_shift = shift.w2step;
    input(cindex,expt,jind).w2min_shift = shift.w2min;
    input(cindex,expt,jind).w2max_shift = shift.w2max;
    input(cindex,expt,jind).shiftstep = shift.step;
    input(cindex,expt,jind).shiftmin = shift.min;
    input(cindex,expt,jind).shiftmax = shift.max;
    input(cindex,expt,jind).shiftb = shiftb;
    
    input(cindex,expt,jind).data_shift_name = data_shift_name;
    input(cindex,expt,jind).pth_data_for_shift = pth_data_for_shift;
    input(cindex,expt,jind).ngr = ngr;
    input(cindex,expt,jind).ni = ni;
    
    input(cindex,expt,jind).thr = thr;
    input(cindex,expt,jind).brem = bneed;
    input(cindex,expt,jind).bins= pulsewb;
    input(cindex,expt,jind).tmini = tmini;
    input(cindex,expt,jind).tmaxi = tmaxi;
    input(cindex,expt,jind).ext= ext;
    input(cindex,expt,jind).wig = wigsb;
    
    input(cindex,expt,jind).pth_wigs = pth_wigs;
    input(cindex,expt,jind).wigsname = wigsname;
    input(cindex,expt,jind).pth_ext = pth_ext;
    input(cindex,expt,jind).extname = extname;
    input(cindex,expt,jind).comment = comment;
    input(cindex,expt,jind).tbac = tbac;
    input(cindex,expt,jind).tfw = tfw;
    input(cindex,expt,jind).split_matin = split_matin;
    input(cindex,expt,jind).reach = reach;
    input(cindex,expt,jind).combine_exposures_FLIMage = combine_exposures;
    input(cindex,expt,jind).tsm = tsm;
end
end


%                         input(cindex,expt,jind).r1s = r1s; input(cindex,expt,jind).r2s = r2s;
%                         input(cindex,expt,jind).r1l = r1l; input(cindex,expt,jind).r2l = r2l;