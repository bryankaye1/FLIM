xmin = 0;
xmax = 3;
xsteps = 1000;
li = 1;

x = xmin:(xmax-xmin)/xsteps:xmax;
y = exp(-x/li);
figure(1);
plot(x,y);

xshift = [0,1];
x1 = 0:.01:3; y1 =  exp(-x1/li);
x2 = 1:.01:3; y2 = exp(-(x2-1)/li);
x3 = 2:.01:3; y3 = exp(-(x3-2)/li);
figure(2); clf; hold on;
plot(x1,y1);
plot(x2,y2);
plot(x3,y3);