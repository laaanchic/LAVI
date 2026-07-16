function [freq_span,ind_above,ind_below] = getFreqSpan(col, i)
ind_below = 0;
ind_above = 0;
% Indices above
flag = 1;
while flag
    tmpi = i-(ind_below+1);
    if tmpi>0 
        if col(tmpi)>0, ind_below = ind_below+1;
        else, flag = 0; 
        end
    else, flag = 0; 
    end
end
% indices below
flag = 1;
while flag
    tmpi = i+(ind_above+1);
    if tmpi<=length(col)
        if col(tmpi)>0, ind_above = ind_above+1;
            else, flag = 0; 
        end
    else, flag = 0; 
    end
end
freq_span = ind_above+ind_below;
