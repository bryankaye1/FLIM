function [w1min, w1max, error] = param(w1est, thr, w1min, w1max, w1step)

    if w1min-w1max~=0
    
 
  
    
    w1estm = w1est<thr;
       
    w1mini = find(w1estm==0, 1, 'first');
    w1maxi = find(w1estm==0, 1, 'last');
    
    w1min = (w1mini-1)*w1step+w1min;
    w1max = (w1maxi+1)*w1step+w1min;
    
        if w1maxi >= length(w1estm) || w1mini <= 1 
            error = 'max or min at edge of marginalized vector';
            break
            
        else

                if w1estm(w1mini-1)> w1est(w1mini) || w1estm(w1maxi+1) > w1estm(w1maxi)
                    error ='not monatonic about max';

                else
                    error ='none';
                end

        end

    else

       error ='none';
    end


    
    
    
    
    
    