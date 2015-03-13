
cyclesmax = 15;
for cycle = 0:cyclesmax-1;
    sn = 3 - floor((cycle)/5);
    str(cycle+1).name = strcat('m',num2str(sn),'-spot',num2str(mod(cycle,5)+1),'.sdt');
    fprintf('%s\n',str(cycle+1).name)
end