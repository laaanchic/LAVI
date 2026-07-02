function prev = dispRMVprev (str,prev)
fprintf(repmat('\b',1,prev))
prev = fprintf(str);
