function [t_before, t_after, pks, trgh] = get_burst_durations(i,j,x,pmtr,ii,filtd,n,pwrr)

pks=[];
trgh=[];
if filtd(ii(x),j(x)) > filtd(ii(x),j(x)+1) && filtd(ii(x),j(x)) > filtd(ii(x),j(x)-1)
    pks = [pks, j(x)]; end
if filtd(ii(x),j(x)) < filtd(ii(x),j(x)+1) && filtd(ii(x),j(x)) < filtd(ii(x),j(x)-1)
    trgh = [trgh, j(x)]; end

go_back = 1; % these short while loops are much faster than "find"
t_before = j(x)*pmtr.downsmpFact-1;

while go_back
    if filtd(ii(x),t_before) > filtd(ii(x),t_before+1) && filtd(ii(x),t_before) > filtd(ii(x),t_before-1)
        pks = [pks, t_before]; end
    if filtd(ii(x),t_before) < filtd(ii(x),t_before+1) && filtd(ii(x),t_before) < filtd(ii(x),t_before-1)
        trgh = [trgh, t_before]; end
    if ~pwrr(i(x),t_before) || t_before==1
        go_back = 0;
    else
        t_before = t_before-1;
    end
end
t_before = t_before + 1;

go_fwd = 1;
t_after = j(x)*pmtr.downsmpFact+1;
while go_fwd
    if filtd(ii(x),t_after) > filtd(ii(x),t_after+1) && filtd(ii(x),t_after) > filtd(ii(x),t_after-1)
        pks = [pks, t_after]; end
    if filtd(ii(x),t_after) < filtd(ii(x),t_after+1) && filtd(ii(x),t_after) < filtd(ii(x),t_after-1)
        trgh = [trgh, t_after]; end
    if ~pwrr(i(x),t_after) || t_after==n
        go_fwd = 0;
    else
        t_after = t_after + 1;
    end
end
t_after = t_after - 1;
pks = sort(pks);
trgh = sort(trgh);

