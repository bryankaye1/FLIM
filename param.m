function [w1minout, w1maxout, error] = param(w1est, thr, w1min, w1max, w1step, sl, sr)

if w1min-w1max~=0
    
    w1estm = w1est<thr;
    
    w1mini = find(w1estm==0, 1, 'first');
    w1maxi = find(w1estm==0, 1, 'last');
    
    w1minout = (w1mini-sl)*w1step+w1min;
    w1maxout = (w1maxi+sr)*w1step+w1min;
    
    error = 'none';
    
    %Make sure search space does not grow%
    if w1min > w1minout
        w1minout = w1min;
    end
    
    if w1max < w1maxout
        w1maxout = w1max;
    end
    
    %Set errors in necessary%
    %     if w1estm(w1mini-1)> w1est(w1mini) || w1estm(w1maxi+1) > w1estm(w1maxi)
    %         error ='not monatonic about max';
    %     else
    %         error ='none';
    %     end
    
    if w1maxi >= length(w1estm)
        error = 1;%'max at edge of marginalized vector';
    end
    
    if w1mini <= 1
        error = 0;%'min at edge of marginalized vector';
    end
    
else
    
    w1minout = w1min;
    w1maxout = w1max;
    error ='0';
    
end







