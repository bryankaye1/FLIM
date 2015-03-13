tempf = load('control_shortlife');
edges = 0:4097;

dlong = tempf.dshort;
hislong = histc(dlong, edges);