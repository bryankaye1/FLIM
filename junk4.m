tempf = load('control_longlife');
edges = 0:4097;

dlong = tempf.dlong;
hislong = histc(dlong, edges);