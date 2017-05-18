% Enter "type" = 1 to read out posterior of lifetime
% Enter "type" = 0 to read matouts to plot # of short lifetimes photons vs longlifetime photons
% enter "int" = low/med/hi for auto plots
close all;
clear;
imin_vec = 30666;%[10090,10106,10122,10138,10154,10170];%10000;%
imax_vec = 30667;%[10105,10121,10137,10153,10169,10185];%10005;%

acquisition_time = [100,20,50,100,20,50];
w2verb = 1;

for fov_ind = 1:length(imax_vec)
    i = imin_vec(fov_ind);
    imax = imax_vec(fov_ind);
    ind =0;
    w1vec = [];
    w1MAPvec = [];
    w2vec = [];
    w2errvec = [];
while i<imax+1
    ind = ind+1;
    [inputm,outputm,flag] = load_mat_data(i,'local',1);
    %If matout/in exists, execute code
    if ~flag
        %%Plot/Print results from short lifetime fitting
        if outputm(1,1,1).w2min~=outputm(1,1,1).w2max
            if inputm(1,1,1).ngr==1
               if w2verb
                   [model,pn,res,sumres] = fitFLIM(outputm,1,['matnum',num2str(i)]);
                   figure(ind); clf;
                   positionVector1 = [.15,.7,.8,.2]; %[left, bottom, width, height]
                   subplot('Position',positionVector1);
                   plot(outputm(1,1,1).w2estx,outputm(1,1,1).w2est);
                   title([strrep(outputm(1,1,1).dataname,'_',' '),...
                       ' (matout',num2str(i),')']);
                   xlabel('lifetime (ns)');
                   ylabel('P(\tau_{long})');
                    
                   positionVector2 = [.15,.3,.8,.3];
                   subplot('Position',positionVector2);
                   plot(pn, '.','MarkerSize',5, 'Color', 'b'); hold on; plot(model,'r');
                   ylabel('Normalized Counts');
                   legend('data','model','Location','NorthEast');
                   yl = ylim;
                   text(1500,yl(2)*.75,['\Sigmares = ',num2str(sumres,'%10.1e')]);
                   
                   positionVector3 = [.15,0.1,.8,.1];
                   subplot('Position',positionVector3);
                   plot(res); hold on;
                   plot(zeros(1,3000),'r--');
                   xlabel('bin number');
                   ylabel('Residue');
                   
                   fprintf('w2Best is %s\n',num2str(outputm(1,1,1).w2Best));
               end
                w2vec = [w2vec outputm(1,1,1).w2Best];
                w2errvec = [w2errvec std_dist(outputm(1,1,1).w2est,outputm(1,1,1).w2estx)];
                
            else
                fprintf('\n matin %s\n', num2str(i));
                for n_index=1:inputm(1,1,1).ngr
                    %fprintf('w2Best is %s\n',num2str(output(1,1,n_index).w2Best));
                    w2vec(ind,n_index) = outputm(1,1,n_index).w2Best;
                end
                delta_t = acquisition_time / length(outputm(1,1,1).ni);
                intensity(ind,:) = outputm(1,1,1).ni./delta_t;
            end
        end
        
        %%Plot/Print results from short lifetime fitting
        if outputm(1,1,1).w1min~=outputm(1,1,1).w1max
            if inputm(1,1,1).split_matin==1
                if outputm(1,1,1).w01MAP ~=0
                w1MAPvec = [w1MAPvec outputm(1,1,1).w1MAP];
                end
                
                if outputm(1,1,1).prBest ~=0
                    figure; subplot(2,1,1);
                    plot(outputm(1,1,1).w1estx,outputm(1,1,1).w1est);
                    reform_dname = strrep(outputm(1,1,1).dataname,'_',' ');
                    title([reform_dname,'. Matout ',num2str(i)]);
                    ylabel('Prob(w1)'); xlabel('\tau (ns)');
                    w1vec = [w1vec outputm(1,1,1).w1Best];
                    w1post_std = std_dist(outputm(1,1,1).w1est,outputm(1,1,1).w1estx');
                    
                    subplot(2,1,2);
                    plot(outputm(1,1,1).prestx,outputm(1,1,1).prest);
                    ylabel('Prob(w01)'); xlabel('Photon Fraction');  
                    fprintf('\n%s. matout %1.0f: w1Best is %1.2f +/-%1.2f',...
                        outputm(1,1,1).dataname,i,outputm(1,1,1).w1Best,w1post_std);
                    
                else
                    fprintf('\n%s. matout %s: FF is 0',outputm(1,1,1).dataname,num2str(i));
                end
            else
                w1cell{ind} = [];
                w01cell{ind}=[];
                w1err_cell{ind} = [];
                w01err_cell{ind} = [];
                for j = 0:outputm(1,1,1).split_matin-1
                    [inputm,outputm,flag] = load_mat_data(i+j);
                    if outputm(1,1,1).prBest ~=0             
                        w1cell{ind} = [w1cell{ind} outputm(1,1,1).w1Best];
                        w01cell{ind} = [w01cell{ind} outputm(1,1,1).prBest];
                        
                        w1err = std_dist(outputm(1,1,1).w1est',outputm(1,1,1).w1estx);
                        w1err_cell{ind} = [w1err_cell{ind} w1err];
                        w01err = std_dist(outputm(1,1,1).prest',outputm(1,1,1).prestx);
                        w01err_cell{ind} = [w01err_cell{ind} w01err];       
                    else
                        w1cell{ind} = [w1cell{ind} 0];
                        w1err_cell{ind} = [w1err_cell{ind} 0];
                        w01cell{ind} = [w01cell{ind} 0];
                        w01err_cell{ind} = [w01err_cell{ind} 0];
                    end
                end
                dataname = strrep(strrep(outputm(1,1,1).dataname,'_',' '),'.sdt','');
                fprintf('%s %1.0f groups: w1= %1.2f +/- %1.2f\n',dataname,...
                    outputm(1,1,1).split_matin,w1cell{ind}(end),w1err_cell{ind}(end));
                delta_t = acquisition_time(ind) / length(outputm(1,1,1).ni);
                intensity(ind,:) = outputm(1,1,1).ni./delta_t;
                         
                figure; title(dataname);
                xlabel('Intensity (cps)');
                yyaxis left
                errorbar(intensity(ind,:),w1cell{ind},w1err_cell{ind},'o');
                ylabel('lifetime (ns)'); 
                ylim([outputm(1,1,1).w1min,outputm(1,1,1).w1max]);
                yyaxis right
                errorbar(intensity(ind,:),w01cell{ind},w01err_cell{ind},'.');
                ylabel('Photon fraction');
                ylim([0,inf]);
                drawnow;
                i = i+j;
            end
        end
    else
        fprintf('DNE\n'); 
        contin = input('DNE continue?');
    end
    i = i+1;
end
end
fprintf('\n');
%%
if inputm(1,1,1).ngr~=1 && inputm(1,1,1).w2min~=inputm(1,1,1).w2max
    xmax = max(max(intensity));
    ymax = max(max(w2vec));
    figure;  hold on;
    %axis([0 xmax 0 ymax]);
    for i = 1:ind
        plot(intensity(i,:),w2vec(i,:),':o');
    end
    title('\tau_{long} sensitiy to intensity, 16 groups per dataset');
    xlabel('Intensity (cps)');
    ylabel('MAP for \tau_{long}');
    %legend('10 sec 2e5','10 sec 2e6','100sec 2e5','100sec 2e6','30sec 2e6','Location','NorthEast');
    legend('2e5 FOV ave','2e5 FOV ave','2e5 FOV ave','1e6 FOV ave',...
        '1e6 FOV ave','1e6 FOV ave','Location','NorthEast');
end

if inputm(1,1,1).ngr==1 && inputm(1,1,1).w1min~=inputm(1,1,1).w1max
    fprintf('\nw1 is %2.2f +/- %2.2f\n',mean(w1vec),std(w1vec));
    fprintf('\nw1MAP is %2.2f +/- %2.2f\n',mean(w1MAPvec),std(w1MAPvec));
elseif inputm(1,1,1).ngr==1 && inputm(1,1,1).w2min~=inputm(1,1,1).w2max 
     fprintf('\nw2 is %2.2f +/- %2.2f\n',mean(w2vec),std(w2vec));
     figure;clf;
     errorbar(1:length(w2vec),w2vec,w2errvec,'.'); hold on;
     plot(1:length(w2vec),ones(length(w2vec),1)*mean(w2vec),'--');
     ylabel('MAP for \tau_{long}');
    % xlabel('25 sec intervals between datapoints');
    % title('Changes in \tau_{long} during time series - donor only');
end

%%
if inputm(1,1,1).ngr~=1 && inputm(1,1,1).w1min~=inputm(1,1,1).w1max
    for j = 1:length(w1cell)
   w1_hi(j) = w1cell{j}(end);
   w1_2hi(j) = w1cell{j}(end-1);
    end
fprintf('short-lifetime of highest intensity is %1.2f +/- %1.2f ns\n',mean(w1_hi),std(w1_hi));
fprintf('short-lifetime of 2nd highest intensity is %1.2f +/- %1.2f ns\n',mean(w1_2hi),std(w1_2hi));   
end
