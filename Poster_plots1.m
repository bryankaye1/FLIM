tempf=load('irf-2013-4-18.mat','-mat','irf', 'bneed', 'pulsewb','wig','tmini','tmaxi');
irf=tempf(1).irf;        %Probability of delay (pdf_IRF)
brem = tempf(1).bneed;   %Number of bins to remove from the 12.5ns to match BH data
bins = tempf(1).pulsewb; %Number of bins that make up 12.5ns
wig = (tempf(1).wig)';
tmini = tempf(1).tmini;
tmaxi = tempf(1).tmaxi;
binskeep = bins-brem;

T = 10;

w2 = 3.77;
ga = irf;
s = T/bins:T/bins:T;

f2 = exp(-s/w2);
f2 = [f2 f2];

f2con = conv(f2,ga);
f2bar = f2con(bins+1:2*bins); 
f2h = f2bar(1:binskeep); %Keep only the appropriate bins
f2h = wig.*f2h;
f2h = f2h/sum(f2h);




T= 12.5;
w00 = 0;
w01 = 0;
w1 = .4;
nps = 10000;
[w00out w01out w02out npho p] = Datarealirf_v2_wig(w00, w01, w1, w2, nps);

a = length(f2h);
xp = 0:10/a:10-10/a;

figure(1); clf; p1 = plot(xp,f2h); xlabel('Time (ns)', 'fontsize', 18); 
ylabel('Probability','fontsize', 18);
title({'Probability of a Photon Arriving at Time t';''}, 'fontsize', 18);
axis([0 10 0 0.001]);
set(p1,'Color','black','LineWidth',1.5);
set(gca,'FontSize',11);

 
 figure(2); clf; p2= plot(xp,p); xlabel('Time (ns)'); ylabel('Photon Counts');
 title('Measured Photons (data)');
 
xlabel('Time (ns)', 'fontsize', 18); 
ylabel('Probability','fontsize', 18);
title({'Probability of a Photon Arriving at Time t';''}, 'fontsize', 18);
%axis([0 10 0 25]);
set(p2,'Color','red','LineWidth',1.5);
set(gca,'FontSize',11);
 
 