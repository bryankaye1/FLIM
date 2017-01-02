clear;

for exday = 0:2
    iminv = [3169,3199,3229];
    imaxv = [3198,3228,3258];
imin = iminv(exday+1); imax = imaxv(exday+1);

ind = 0;
for i = imin:1:imax
    
    clear output intb ffP cintpr
    try
        nstr = strcat('matout',num2str(i),'.mat');
        tempf=load(nstr,'-mat','output','prBestmat','w1Bestmat','w2Bestmat','eltime');
    catch exception
        try
            nstr = strcat('Y:\Users\bkaye\cluster\matout\matout',num2str(i),'.mat');
            tempf=load(nstr,'-mat','output','prBestmat','w1Bestmat','w2Bestmat','eltime');
        catch
            continue
        end
    end
    ind =ind+1;
    output = tempf.output(:,:,:);
    dataname = output.dataname;
    %fprintf('pth is %s\n,',output(1,1,1).pth_data);
    % fprintf('%s\n',dataname);
    x = output(1,1,1).ni;%./output(1,1,1).sinti;
    al = output(1,1,1).w1Best./output(1,1,1).w2Best;
    ep = output(1,1,1).sinti;
    for j = 1:length(output)
        y(j) = output(1,1,j).prBest;
        cintpr(j) = findci(output(1,1,j).prest,output(1,1,j).prestx,.682)/2;
        xerr = sqrt(x);
    end
    
    if min(x)==0
        istart = find(x==0,1,'last')+1;
        x = x(istart:end);
        y = y(istart:end);
        cintpr = cintpr(istart:end);
        xerr = xerr(istart:end);
        ep = ep(istart:end);
    end
    
    
    sumres2 =inf;
    bmin = 0;
    bmax = mean(x);
    amax = max(y)*4;
    amin = amax*.1;
    if amin == amax || bmin ==bmax
        modelb = y;
        ab = pi;
        bmonb = pi;
    else
        for a = amin:(amax-amin)/100:amax
            a3 = 1+a*(al-1);
            for bmon = bmin:(bmax-bmin)/100:bmax;
                model = (a*al/a3)*(1-(ep*bmon)./x);
                res = (y-model).*(x.^1);
                sres = sum(res.^2);
                if sres < sumres2
                    sumres2 = sres;
                    ab = a;
                    bmonb = bmon;
                    modelb = model;
                end
            end
        end
        if bmonb == bmin || bmonb==bmax
            fprintf('bmonb is %f bmin/max is %f\n',bmonb,bmin);
        end  
        if ab == amin || ab==amax
            fprintf('ab is %f amin/max is %f\n',ab,amin);
        end
    end
%     figure(ind); clf; plot(x,y,'-o'); hold all; plot(x,modelb, 'r');
%     axis([0 max(x)*1.1 min(y)*.8 max(y)*1.2]);%plot(x,modelb1, 'g');
%     errorbar(x,y,cintpr,':','Color',[.7 .7 .9]);
%     xlabel('intensity'); ylabel('Fret Fraction');
%     ti = strcat(num2str(dataname(1:end)),'  a=',num2str(ab,3),'  bmon=',num2str(32*bmonb,2));
%     title(ti);
%     fprintf('%s\n',ti);

    bpol(ind) = sum(x.*y./(ep*a*al));
    tub(ind) = sum((x./ep));
    pfmeas(ind) = ab;
end
ave=1;
%%
if ave ==0
    pfm = pfmeas;%pfmeas(end-5:end);%-pfmeas(1);
    sresB= inf;
    for ex=10:30
        acc = 300;
        acclab = acc*.66;
        donor = 250;
        %ex = 15;
        acv = .5:.5:3;
        j = 0; clear m;
        for i = 1:6
            j = j+1;
            m(j) = acclab*acv(j) / (acc*acv(j)+1*donor+40*ex);
        end
        m = [0,m];
        pfpred = 2*m;
        res = ((pfm-pfpred)./sqrt(pfpred)).^2;
        sres = sum(res);
        if sres < sresB
            pfpredB = pfpred;
            exB=ex;
            sresb=sres;
        end
    end
    figure(ind+1); clf; hold all;
    plot(pfm,pfpredB,'-o');
    plot(pfpredB,pfpredB,'r');
    axis([0 max([pfm,pfpred])*1.1 0 max([pfm,pfpred])*1.2]);
    xlabel('predicted pFRET'); ylabel('measured pFRET');
    
    figure(ind+2); clf; hold all;
    %plot(length(pfm),pfm,'-o');
    f = fit((1:length(pfm))',pfm','poly1');
    plot(f,1:length(pfm),pfm);
    axis([0 length(pfm) 0 max(pfm)*1.2]);
end


if ave==1
    bpolr = reshape(bpol,6,5);
    bperr=std(bpolr/max(mean(bpolr,1)),1)/sqrt(6);
    bpy = mean(bpolr,1)/max(mean(bpolr,1));
    bpx = 2.^(1-length(bpy):0);
    figure(ind+1+exday*3); clf; hold all;
    f = fit(bpx',bpy','poly1');
    plot(f,bpx,bpy);
    errorbar(bpx,bpy,bperr,'.','Color',[.7 .7 .9]);
    title('Pol by FRET');
    
    tubr = reshape(tub,6,5);
    tuberr=std(tubr/max(mean(tubr,1)),1)/sqrt(6);
    tuby = mean(tubr,1)/max(mean(tubr,1));
    tubx = 2.^(1-length(tuby):0);
    figure(ind+2+exday*3); clf; hold all;
    f = fit(tubx',tuby','poly1');
    plot(f,tubx,tuby);
    errorbar(tubx,tuby,tuberr,'.','Color',[.7 .7 .9]);
    title('Pol by int');
    
    figure(ind+3+exday*3); clf; hold all;
    f = fit(tuby',bpy','poly1');
    plot(f,tuby,bpy);
    errorbarxy(tuby,bpy,tuberr,bperr,{'.', 'b', 'b'});
    title('Pol from FRET vs Pol from int');
end



old =0;
if old ==1
    
    bpolr = reshape(bpol,3,5);
    mbpol = mean(bpolr,1);
    bpolmin = max(mbpol)/2^(length(mbpol)-1);
    intx = bpolmin*2.^(0:length(mbpol)-1);
    figure(ind+1); clf; hold all;
    plot(intx,mbpol,'bo'); plot(intx,intx,'r'); %errorbar(intx,intm,intv);
    title('Pol by FRET');
    
    tubr = reshape(tub,3,5);
    mtub = mean(tubr,1);
    tubmin = max(mtub)/2^(length(mtub)-1);
    tubx = tubmin*2.^(0:length(mtub)-1);
    
    figure(ind+2); clf; hold all;
    plot(tubx,mtub,'bo'); plot(tubx,tubx,'r'); %errorbar(intx,intm,intv);
    title('polymer by intensity');
    
    figure(ind+3); clf; hold all; plot(tubx,mbpol,'bo'); plot(tubx,tubx,'r');
    title('polymer by FRET vs Polymer by intensity');
end
end