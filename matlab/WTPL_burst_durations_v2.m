function [t_before, t_after] = WTPL_burst_durations_v2(f_ind,B, wtpl)
% wtpl is a binary vector 
x0 = B(2);
go_back = 1;
t_before = x0-1;
while go_back
    if ~any(wtpl(f_ind,t_before)) || t_before==1
        go_back = 0;
    else
        t_before = t_before-1;
    end
end
t_before = t_before + 1;
    
go_fwd = 1;
t_after = x0+1;
while go_fwd
    if ~any(wtpl(f_ind,t_after)) || t_after==length(wtpl)
        go_fwd = 0;
    else
        t_after = t_after + 1;
    end
end
t_after = t_after - 1;