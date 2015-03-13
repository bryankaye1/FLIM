function [errPR, errw02, errw2, errw1] = errorcheck(prBesti,prest,prBest,w02Besti,w02est,w02Best,w2Besti,w2est, w1Besti,w1est,l1)

errPR = 0;
errw02 = 0;
errw1 = 0;
errw2 = 0;

if prBesti == length(prest) || prBesti == 1
    if prBest ~= 0 && prBest ~=1 && range(prest)~=0
        errPR = 1;
    end
end

if w02Besti == length(w02est) ||w02Besti == 1
    if w02Best ~= 0 && w02Best ~=1 && range(w02est)~=0
        errw02 = 1;
    end
end

if w2Besti == length(w2est) || w2Besti == 1
    if range(w2est)~=0
        errw2 = 1;
    end
end

if w1Besti == length(w1est) || w1Besti == 1
    if range(w1est)~=0
        errw1 = 1;
    end
end


end