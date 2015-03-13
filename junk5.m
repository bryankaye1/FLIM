%tempfs = load('control_shortlife');
edges = 0:4095;
% 
% dshort = tempfs.dlong; 
hisshort = histc(dshort, edges);

% tempfl = load('control_longlife');
% edges = 0:4097;
% 
% save('control_longlife', 'hislong');  
% 
% 
% 
% dlong = tempfl.dlong;
 hislong = histc(dlong, edges);
% 
plot(hisshort,'r'); hold on; plot(hislong,'b');
legend('short','long');
% 
% save('short_his', 'hisshort');  
% 
% save('long_his', 'hislong');  