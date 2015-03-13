%% Clears and Comments
clear;
%clc; 
clf;

%Inputs: Donor lifetime. 
%Outputs: Background, FRET lifetime, FRET and Non-FRET fractions
%V8: added w1min

%Need to add:
%1) For window factor, don't use best in each dimension, but instead use best in parameter space

%2) Marginalize to find frac_est, then apply scale factor

%% Load IRF
tic

tempf=load('irf-2013-4-18.mat','-mat','irf', 'bneed', 'pulsewb','wig');
irf=tempf(1).irf;        %Probability of delay (pdf_IRF)
brem = tempf(1).bneed;   %Number of bins to remove from the 12.5ns to match BH data
bins = tempf(1).pulsewb; %Number of bins that make up 12.5ns
wig = (tempf(1).wig)';
binskeep = bins-brem;


%% Initial Conditions
%%% USED FOR marginalization at end
%Insert hist data to get p, search over w0 and w1
w1step = .001; w1min=0.30; w1max = .58; %w1min must be an integer multiple of w1step
fracstep = 0.0005;
w00step = fracstep; w00min=0.005; w00max = .015;
w02step = fracstep; w02min = 0.90; w02max = 1;

%%%USED FOR BAR GRAPHS
% w1step = .01; w1min=0.30; w1max = .58; %w1min must be an integer multiple of w1step
% fracstep = 0.001;
% w00step = fracstep; w00min=0.005; w00max = .015;
% w02step = fracstep; w02min = 0.90; w02max = 1;

T = 12.5;
w00 = 0.01;
w01 = .05;
w1 = .4;
w2 = 3.822;
nps = 1000000;
[w00out w01out w02out npho p] = Datarealirf_v2_wig(w00, w01, w1, w2, nps);
%[p file_name] = read_data_14ns;

w2 = 3.822;

%% Calculate Post 

ga = irf;
s = T/bins:T/bins:T;

loglike = -pi*ones(round(1+(w1max-w1min)/w1step),round(1+(w02max-w02min)/w02step),round(1+(w00max-w00min)/w00step));
  sizell = size(loglike,1)*size(loglike,2)*size(loglike,3)/1000; 
  time_est= ceil(.1294*sizell);
  fprintf('\n\n%1.0fk gridpoints. Runtime %1.0f seconds\n',sizell,time_est);
    
f2 = exp(-s/w2);
f2 = [f2 f2];

f2con = conv(f2,ga);
f2bar = f2con(bins+1:2*bins); 
f2h = f2bar(1:binskeep); %Keep only the appropriate bins
f2h = f2h/sum(f2h);
    
for w1i = w1min:w1step:w1max
  %w1i     
clear f;   
f = exp(-s/w1i);
f = [f f];
fcon = conv(f,ga);
fbar = fcon(bins+1:2*bins);
fh = fbar(1:binskeep);
fh = fh/sum(fh);

for w02i = w02min:w02step:w02max
for w00i = w00min:w00step:w00max %Now we just have to search the space
 
    if w00i + w02i > 1
    loglike(round(1+(w1i-w1min)/w1step),round(1+(w02i-w02min)/w02step),round(1+(w00i-w00min)/w00step)) = -inf;
    else
    loglike(round(1+(w1i-w1min)/w1step),round(1+(w02i-w02min)/w02step),round(1+(w00i-w00min)/w00step)) =...
        sum(log(wig.*(w00i/binskeep+(1-w00i-w02i)*fh + w02i*f2h)).*p); %
    end          

end
end
end
%There are NANs in the loglike whenever a bin has w00i = w02i = 0 and fh =0
%for a certain bin, which is gives log(0)= -inf. If this is then multiplied
%by p > 0, this gives a NAN (0*-inf =nan).
b = loglike(find(loglike==-pi));
if b~=0
    disp('You are fucked, there is a pi in there');
end

loglike2 = loglike - max(max(max(loglike))); %If there are imaginary components- BAD
like = exp(loglike2);
%%%%%% Priors %%%%%%

pri = zeros(round(1+(w1max-w1min)/w1step),round(1+(w02max-w02min)/w02step),round(1 + (w00max-w00min)/w00step)); %Ini Properly
pri(:,1,:)= 1;
post = like.*pri; %Change all likes below to post

%% Marg and Plots

    %%%%%% Marginalize
w00est = squeeze(sum(sum(like,1),2));  %Can replace with trapz's later. 
w00estx = w00min:w00step:w00max;            %See last part of desc. section of trapz for problem
[c i]= max(w00est);
w00Best = (i-1)*w00step+w00min;

w1est = squeeze(sum(sum(like,2),3));
w1est= w1est/sum(w1est);
w1estx = w1min:w1step:w1max;
[c i]= max(w1est);
w1Best = (i-1)*w1step+w1min;

w02est = squeeze(sum(sum(like,1),3));
w02estx = w02min:w02step:w02max;
[c i]= max(w02est);
w02Best = (i-1)*w02step+w02min;

            %%%%%%%%%%%%%%%% PLOTTING %%%%%%%%%%%%%
 %% 3D surf Plot           
%%%Plot in full window (no plot assitant) on small screen
 [X,Y] = meshgrid(w02min:w02step:w02max,w1min:w1step:w1max); 
Z = squeeze(sum(like,3))/sum(sum(squeeze(sum(like,3))));
figure(2);clf;
bar4viacolor(struct('Z',Z));
xlabel('FRET Fraction','fontsize',18);ylabel('FRET Lifetime','fontsize',18);
zlabel('Probability of Parameters','fontsize',18);title('1M Photons','fontsize',18);


%set(gca, 'XTick', w02min:w02step:w02max);
xvec = w02min+(w02max-w02min)/10:(w02max-w02min)/10:w02max;
xvec = floor(100*xvec)/100;
set(gca,'XTickLabel',xvec,'fontsize',14);
xlabh = get(gca,'XLabel'); xpos = get(xlabh,'Position') + [20 0 0];
set(xlabh,'Position',xpos);


yvec = w1min:(w1max-w1min)/3:w1max; yvec(end) = .6;
yvec = floor(100*yvec)/100;
set(gca,'YTickLabel',yvec, 'fontsize',14);
ylabh = get(gca,'YLabel');
set(ylabh,'Position',[-274.577 366.018 0.204]);

zlabh = get(gca,'ZLabel'); zpos = get(zlabh,'Position') + [0 0 0.005]; %Pos[]= Xpos,Ypos,Zpos in units of data
set(zlabh,'Position',zpos);

%titleh = get(gca, 'Title'); set(titleh, 'Position', [-232.88625885254652 366.1798006578003 0.2388925311667347]);

%%Move and shrink color bar%%
h_bar = findobj(gcf,'Tag','Colorbar');
initpos = get(h_bar,'Position');
initfontsize = get(h_bar,'FontSize');

set(h_bar,'Position',[initpos(1)*1.08 initpos(2)*1.5 initpos(3) initpos(4)*.7]);%Pos[] = [xpos ypos width hieght]
   %...'FontSize',initfontsize) to make fontsize of color bar different


%% Marginalized Plot
% figure(1); clf; surf(X,Y,Z); xlabel('N_{0}');ylabel('\tau(ns)');colormap(jet);
% zlabel('Probability(N_{0},\tau)');title('Simulated Data: 100M Photons');
w02est = w02est/sum(w02est);
figure(3); clf; p3=plot(w02estx,w02est);
ylabel('Probability of FRET Fraction','fontsize', 20); xlabel('FRET Fraction','fontsize', 20);
title('1M Photons','fontsize', 20);

set(gca,'FontSize',12);
set(p3,'Color','black','LineWidth',2);
%figure(2); clf; plot(w1estx,w1est); xlabel('\tau_{F}(ns)'); ylabel('Prob(\tau_{F})'); title('Probability(\tau_{F}|data)');
%figure(1); clf; plot(w00estx,w00est); xlabel('w00'); ylabel('P(w00)');
%bar4viacolor(Z);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
             
%% Pop estimates
           
%%%%% FOR HIGH PHOTON NUMBER, this is equal to the relative fraction
w01Best = 1 - w02Best - w00Best;
w02 = 1 - w01 - w00;
%Psig = w01/(w01 + w00);
%Psigout = w01out/(w01out+w00out);
%Psest = w01Best/ (w01Best+w00Best);

%%%%% Find "a/b" factor %%%%%

%%%%% find b %%%%%%
clear f2;
f2 = exp(-s/w2);
f2 = [f2 f2];

f2con = conv(f2,ga);
f2bar = f2con(bins+1:2*bins); 
f2h = f2bar(1:binskeep); s2t = sum(f2bar); %Keep only the appropraite bins
s2 = sum(f2h);
b = s2/s2t;

%%%%%% find a %%%%%%
clear f;
f = exp(-s/w1Best);
f = [f f];
fcon = conv(f,ga);
fbar = fcon(bins+1:2*bins); s1t = sum(fbar);
fh = fbar(1:binskeep);
s1 = sum(fh);
a = s1/s1t;

r = (w01Best/w02Best)*(b/a);
rpop = r*w2/w1Best;
fracestr = 100*rpop/(rpop+1); %Def works, when testing, makes sure grid spacing is sufficiently small %see v7 for testing

%% Print outs

% %fprintf('\n There were %3.0f photons in the data set \n', npho);
% fprintf(' w00Best %3.4f  \n', w00Best);
% fprintf(' w02Best %3.4f  \n', w02Best);
fprintf(' w1Best%3.4f   w1 %3.4f \n', w1Best, w1);
fprintf(' w01Best %3.4f   w01out %3.4f \n', w01Best, w01out);
%fprintf(' w1est  %3.4f     w1act %3.4f\n', w1Best, w1);



