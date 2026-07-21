function txt = showCurrnetTime
c=clock; 
hr = num2str(round(c(4)));
min = num2str(round(c(5)));
if length(min)<2; min = ['0' min]; end
sec = num2str(round(c(6)));
if length(sec)<2; sec = ['0' sec]; end
txt = [hr ':' min ':' sec];
end