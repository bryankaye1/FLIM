%generate FRET and intensity data for non-uniform monomer concentration
clear;

% Initial stuff
addnoise = 1; %set to 1 if you want to add noise to FRET and intensity
ind = 1; % for setting the figure number
pf = .14; %Probability of donor engaged in FRET, given it is in polymer.
ep = 5e2; % intensity normalization
al = 1/3.7; %relative intensity of donors engaged in FRET to donors not engaged

%Set Polymer
Npol_out = [1:.05:2,3:1.5:9]; %polymer at the interface and out of the spindle
Npol_in = [12:.2:13]; %polymer in the spindle
Npol = [Npol_out,Npol_in]; 

%Set Monomer %%%SEBASTIAN LOOK HERE TO PLAY WITH MONOMER CONCENTRATION%%%
mon_amp_in = 0.5; %amplitude of monomer in the spindle
mon_amp_out = 1;
Nmon = 3*[mon_amp_out*ones(1,length(Npol_out)),mon_amp_in*ones(1,length(Npol_in))];

%FRET and intensity given monomer and polymer
fret = pf*Npol./(Npol + Nmon) ;
int = ep * (Nmon + (1-pf)*Npol + al*pf*Npol);

%Add noise
if addnoise
    %Add noise to FRET
    noise_std = (sqrt(int(1))./sqrt(int)) .*0.005; %assume noise scales as 1/sqrt(N)
    fret_noise = randn(1,length(fret)).*noise_std; %gaussian noise with std noise_std
    fret = fret+fret_noise; % add noise to FRET signal.
    %enforce that FRET is between 0 and 1. 
    fret(fret<0)=0;
    fret(fret>1)=1;
    
    %add noise to intensity
    int = poissrnd(int);
end

%Fit data with spatially uniform monomer model
[pfm(ind),nmon(ind),fresult] = fit_pf_nm(int,fret,noise_std,al);

%Plot it
figure(ind); clf; hold on; 
errorbar(int,fret,noise_std,'ko');%plot the FRET vs intensity with errorbars
xmodel = 0:max(int)*.001:max(int)*1.2;%vector for plotting the fit
plot(xmodel,pf_nm_pred(al,pfm(ind),nmon(ind),xmodel),'--','Color',[0.6,0.6,0.6]); %plot the fit
title('simulated data');
xlabel('intensity (au)');
ylabel('FRET fraction');