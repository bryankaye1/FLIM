function [w1MAP,w2MAP,w01MAP,w02MAP] = getMAP(like,w1estx,w2estx,prestx,w02estx)
%Find MAP from likelihood

[w1ind,w2ind, w02ind, w01ind]= ind2sub(size(like),find(like==max(max(max(max(like))))));

stepsize = @(vec) ( vec(end)-vec(1) ) / (length(vec)-1 );
xmax = @(ind,vec) (ind-1) * stepsize(vec) + vec(1);

w1MAP = xmax(w1ind,w1estx);
w2MAP = xmax(w2ind,w2estx);
w01MAP = xmax(w01ind,prestx);
w02MAP = xmax(w02ind,w02estx);

end

%     if w1Best ~=3 && w1Best ~=1
%     fprintf('%1.2f w1 / %1.2f w01 / %1.2f w02 matin%1.0f %s\n', ...
%         w1Best, w01Best, w02Best,i,outm(1,1,1).dataname);
%     end
%
%     else
%         fprintf('multiple MAPS. matin%1.0f\n', i);
%     end

