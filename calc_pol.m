function [Npol_ave, Npol_std] = calc_pol(al,pf,FRET_ave,FRET_std,intensity)

imax = 100000;
for i = 1:length(FRET_ave)
    FRET = normrnd(FRET_ave(i),FRET_std(i),1,imax);
    Npol_sams = FRET.*intensity(i) ./ (pf*(1+(al-1).*FRET));
    Npol_ave(i) = mean(Npol_sams);
    Npol_std(i) = std(Npol_sams);
end

%figure(1); 
%hist(Npol_sams,100);
%pause;


end
