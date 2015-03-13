clear; clc;
%This program shows the output of the IRF conv with signal, then moded, then moving times. 
%The program also show the measured IRF (which has moved times) convolved with signal

%"p" corresponds to the action of moving times

%Need to add: missing time simulation 
%-ei if you put the missing time in the wrong spot


for i = 1:8 %loop to see how difference between two methods differs with grid spacing

xstep = .01*2^(-i+1); %Time point spacing
xmax = 1;
nxs = round(xmax/xstep); %Number of time points
x = xstep: xstep: xmax;

%IRF
irf = zeros(1, nxs);
irf(round(.3*nxs)) = 1;
irf(round(.5*nxs)) = .5;
irf(round(.95*nxs)) = .1;
%irf = (xstep:xstep:xmax); irf = irf.*irf; %squared t IRF
figure(1); plot(irf); title('irf');


%PIRF
pirf = [irf(round(.8*nxs+1): nxs) irf(1:round(.8*nxs))];
figure(2); clf;  plot(pirf); title('pirf');

 
%SIGNAL
sig = exp(-x);
sig = [sig sig];
figure(3); clf; plot(sig); title('sig');


%SIGNAL CONVOLVED WITH IRF
cis = conv(irf,sig);
figure(4); clf; plot(cis); title('sig conv irf');


%MODED -SIGNAL CONVOLED WITH IRF-
mcis = cis((nxs+1):2*nxs);
figure(5); clf; plot(mcis); title('moded sig conv irf');


%MOVE TIMES AROUND OF 'MODED -SIGNAL CONV IRF-'
pmcis = [mcis(round(.8*nxs+1): nxs) mcis(1:round(.8*nxs))];
pmcis = pmcis/sum(pmcis);
figure(6); clf; plot(pmcis); title('P of moded sig conv irf');


%MEASURED IRF (MOVED TIMES IRF) CONVOLVED WITH SIGNAL
csp = conv(pirf,sig);
figure(7); clf; plot(csp); title('sig conv pirf');

%MOD MEASURED IRF CONV SIGNAL
mcsp = csp((nxs+1):2*nxs);
mcsp = mcsp / sum(mcsp);
figure(8); clf; plot(mcsp); title('moded sig conv pirf');

%CALCULATE DIFFERENCE OF TWO SCENARIOS AS A FUNCTION OF TIME SPACING
meanans(i) = mean(abs(mcsp-pmcis)./mcsp);
maxans(i) = max(abs(mcsp-pmcis)./mcsp);
end

figure(9); clf; plot(log10(maxans)); title('Max diff vs bins -Orders of Mag');



 
%  
% xt1 = .1:.1:2;
% sigt1 = exp(-xt1);
%  
% sigt2 = [sig sig];
%  
% ct1 = conv(irf, sigt1);
% ct2 = conv(irf, sigt2);
%  
% t1 = ct1(1:10) + ct1(11:20) + ct1(21:end);
% t2 = ct2(1:10) + ct2(11:20) + [ct2(21:end) ct2(end)];
% t3 = cis(1:10) + cis(10:end);
% t4 = ct2(10:19);
% t5 = ct2(11:20);
%  
% t1= t1/sum(t1);
% t2 = t2/sum(t2);
% t3= t3/sum(t3);
% t5 =t5/sum(t5);
%  
% figure(1); clf; plot(cis);
% figure(2); clf; plot(t2);
% figure(3); clf; plot(t5);
% figure(4); clf; plot(ct2/sum(sig))