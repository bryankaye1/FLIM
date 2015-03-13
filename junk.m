pr = .2;
w02 = .95;
w1 = .5;
w2 = 3.8;
nps = 5000000;

[w00out, w01out, w02out, npho, dt] = SimData_v1(pr, w02, w1, w2, nps);

figure(1); plot(dt);