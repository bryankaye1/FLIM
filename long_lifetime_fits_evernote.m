% Enter "type" = 1 to read out posterior of lifetime
% Enter "type" = 0 to read matouts to plot # of short lifetimes photons vs longlifetime photons
% enter "int" = low/med/hi for auto plots
%close all;
clear;
i = 7556;%10000;%7561%10005 %Enter matstart here
imax = 7560;%10004;%7565%10009 %Enter matend here 
ind =0;
%Enter acquisiton time here. THis is to convert total photons to cps
acquisition_time = [10,10,100,100,30];


while i<imax+1
    ind = ind+1;
    [input,output,flag] = load_mat_data(i);
    %If matout/in exists, execute code
    if ~flag
        %%Plot/Print results from short lifetime fitting
        if output(1,1,1).w2min~=output(1,1,1).w2max
            if input(1,1,1).ngr==1
                [model,pn,res,sumres] = fitFLIM(output,1,['matnum',num2str(i)]);
                figure(ind); clf;
                positionVector1 = [.15,.7,.8,.2]; %[left, bottom, width, height]
                subplot('Position',positionVector1);
                plot(output(1,1,1).w2estx,output(1,1,1).w2est);
                title([strrep(output(1,1,1).dataname,'_',' '),...
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

                fprintf('w2Best is %s\n',num2str(output(1,1,1).w2Best));
                
            else
                fprintf('\n matin %s\n', num2str(i));
                for n_index=1:input(1,1,1).ngr
                    %fprintf('w2Best is %s\n',num2str(output(1,1,n_index).w2Best));
                    w2vec(ind,n_index) = output(1,1,n_index).w2Best;
                end
                delta_t = acquisition_time(ind) / length(output(1,1,1).ni);
                intensity(ind,:) = output(1,1,1).ni./delta_t;
            end
        end
        
        %%Plot/Print results from short lifetime fitting
        if output(1,1,1).w1min~=output(1,1,1).w1max
            if input(1,1,1).split_matin==1
                if output(1,1,1).prBest ~=0
                    figure; plot(output(1,1,1).w1estx,output(1,1,1).w1est); title(['w1 for matout ',num2str(i)]);
                    fprintf('\nmatout %s: w1Best is %s',num2str(i),num2str(output(1,1,1).w1Best));
                    figure; plot(output(1,1,1).prestx,output(1,1,1).prest); title(['w01 for matout ',num2str(i)]);
                end
            else
                fprintf('\n\n matin %s', num2str(i));
                w1vec{ind} = [];
                w01vec{ind}=[];
                for j = 0:output(1,1,1).split_matin-1
                    [input,output,flag] = load_mat_data(i+j);
                    
                    if ~flag
                        if output(1,1,1).prBest ~=0
                            fprintf('\n int_group %s: ', num2str(j+1));
                            fprintf('w1Best is %s, w01Best is %s',...
                                num2str(output(1,1,1).w1Best),num2str(output(1,1,1).prBest));
                            w1vec{ind} = [w1vec{ind} output(1,1,1).w1Best];
                            w01vec{ind} = [w01vec{ind} output(1,1,1).prBest];
                        end
                    else
                        v1vec{ind} = [w1vec 0];
                    end
                end
                figure;
                yyaxis left
                plot(w1vec{ind},'o');
                yyaxis right
                plot(w01vec{ind},'+');
                title(strcat('w1Best for', num2str(i)));
                i = i+j;
            end
        end
    else
        fprintf('DNE\n');
        
    end
    
    %     if output(1,1,1).prmin~=output(1,1,1).prmax
    %
    %     end
    i = i+1;
end
fprintf('\n');
%%
if input(1,1,1).ngr~=1
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
    legend('10 sec 2e5','10 sec 2e6','100sec 2e5','100sec 2e6','30sec 2e6','Location','NorthEast');
    
end



%%
