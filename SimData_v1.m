function [w00out, w01out, w02out, npho, dt] = SimData_v1(pr, w02, w1, w2, nt)

%%Add wiggles vector, deconvolve wiggles from IRF
T=12.58;

%%%%%%%% Load in IRF information %%%%%%%%
                tempf=load('Z:\bkaye\Bayes_2013\mat_files\current.mat','-mat','bneed',...
                    'pulsewb', 'irf', 'irfsim', 'irf_pdf', 'tmini', 'tmaxi', 'ext','irfname','wigsb','gab');

brem = tempf(1).bneed;   %Number of bins to remove from the 12.5ns to match BH data
bins = tempf(1).pulsewb; %Number of bins that make up 12.5ns
tmini = tempf(1).tmini;
tmaxi = tempf(1).tmaxi;
ext = tempf(1).ext;
wig = tempf(1).wigsb;
ga = tempf(1).gab;
binskeep = bins-brem;


[irftimes]=pdf_builder(ga);


%%%%%% Convert from pr and w02 %%%%%%
w01 = pr*w02*w1/(w2*(1-pr));%(w1/w2)* w02*pr;
w00 = 1 - w01 - w02;

%%%%%%%% Simulate Data %%%%%%%%
np0 = round(w00*nt);
np1 = round(w01*nt);
np2 = round((1-w00-w01)*nt);

%Generate Signal Photons
rng('shuffle');
r1 = rand(1,np1);
dts1 = -w1.* log(1 - r1); %dts1 is emission time of signal photons from w1
r2 = rand(1,np2);
dts2 = -w2.* log(1 - r2); %dts2 is emission time of signal photons from w2

%Generate background photons
npb = round(np0); %nbd is number of background photons
dtb = T * rand(1,npb); % dtb is the emission time of background photons

%Generate Proper delay times of each photon
delays1=irftimes(randi(length(irftimes),1,length(dts1)));
delays2=irftimes(randi(length(irftimes),1,length(dts2)));
delayb=irftimes(randi(length(irftimes),1,length(dtb)));

%add delays and mod data
ts1 = mod(dts1+delays1, T);
ts2 = mod(dts2+delays2, T);
tb = mod(dtb+delayb, T);

%Histogram data
edges = T/bins:T/bins:T; %Create equally spaced time bins
s1 = histc(ts1,edges); %Histogram of data
s2 = histc(ts2,edges);
b = histc(tb,edges);

%%Remove time bins from the times at which you cannot detect photons
binskeep = bins - brem;
s1 = s1(1:binskeep); s1 = round(s1.*wig); 
s2 = s2(1:binskeep); s2 = round(s2.*wig);
b = b(1:binskeep); b = round(b.*wig);


dt = s1+s2+b; %add signal and background photons
npho = sum(dt); %Total photons left after we remove undetecable times
w00out = sum(b)/npho; %Fraction of background photons
w01out = sum(s1)/npho; %Fraction of w1 photons
w02out = sum(s2)/npho; %Fraction of w2 photons
end

