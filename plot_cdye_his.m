sh = load('sh-his.mat');
ysh = sh.info.his.his;
figure(1); clf; plot(log(ysh/sum(ysh)),'b'); hold all;

sm = load('sm-his.mat');
ysm = sm.info.his.his;

plot(log(ysm/sum(ysm)),'g');

sl = load('sl-his.mat');
ysl = sl.info.his.his;

plot(log(ysl/sum(ysl)),'r');


lh = load('lh-his.mat');
ylh = lh.info.his.his;
figure(2); clf; plot(log(ylh/sum(ylh)),'b'); hold all;

lm = load('lm-his.mat');
ylm = lm.info.his.his;

plot(log(ylm/sum(ylm)),'g');

ll = load('ll-his.mat');
yll = ll.info.his;

plot(log(yll/sum(ylh)),'r');