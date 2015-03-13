figure(20); clf;  hold all; 
xlabel('FRET FRAC'); ylabel('Fret Lifetime');
title('Taxol-Fresh Extract: IRF828pm (Taxol) Donor-Acceptor Data used for IRF Shift');
axis([0 .12 0 1.4]);
for i = 1:6
y = cell2mat(w1Bestmat(:,2,i));
x = cell2mat(prBestmat(:,2,i))./(1+cell2mat(prBestmat(:,2,i)));
scatter(x,y);
end
legend('Slide1 No Thresh','Slide2 No Thresh',...
    'Slide1 1M Thresh','Slide2 1M Thresh', 'Slide1 2M Thresh','Slide1 2M Thresh',...
    'Location','SouthEast');


figure(21); clf; hold all; 
xlabel('FRET FRAC'); ylabel('Fret Lifetime'); 
title('Taxol-Fresh Extract: IRF235pm (No Taxol)donor only Data used for IRF Shift');
axis([0 .12 0 1.4]);
for i = 1:6
y = cell2mat(w1Bestmat(:,1,i));
x = cell2mat(prBestmat(:,1,i))./(1+cell2mat(prBestmat(:,2,i)));
scatter(x,y);
end
legend('Slide1 No Thresh','Slide2 No Thresh',...
    'Slide1 1M Thresh','Slide2 1M Thresh', 'Slide1 2M Thresh','Slide1 2M Thresh',...
    'Location','SouthEast');